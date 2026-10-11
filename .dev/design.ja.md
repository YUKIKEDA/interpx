# `interpx` Architecture & Design Specification

> **Pure Rust による `scipy.interpolate` 主要補間クラス網羅の補間ライブラリ。**

基準とする SciPy のバージョンは **1.18** で固定する（数値一致テスト・網羅表の定義はすべてこの版に従う）。

---

## 1. ビジョンと設計原則 (Vision & Design Principles)

### 1.1 背景と目的
Rust の科学計算エコシステムにおける最大の課題は「断片化（fragmentation）」である。1次元補間、グリッド補間、散乱点補間が個別のクレートに分かれ、API や型システムが統一されていない。特に不規則散乱データの多次元補間（`LinearNDInterpolator`, `griddata`）は実質的な空白地帯となっている。

`interpx` は、**`scipy.interpolate` の主要補間クラス（§5 の網羅表）を網羅するオールインワンのデファクトスタンダード**を目指す。網羅の範囲外は §5.3 の対象外リストに明記し、ロードマップ（§9）で扱う。

### 1.2 コア設計原則
1. **Pure Rust (C/Fortran 依存ゼロ)**: `cargo build` のみでどこでもビルド可能（クロスコンパイル対応）。
2. **静的型安全性**: `const D: usize` によるコンパイル時次元最適化とスタック割り当てを優先。
3. **Rust イディオムと高スループットの両立**: `Option` による安全な領域外判定、メモリ再利用（ゼロアロケーション）、`rayon` による並列バッチ評価。
4. **モジュラー**: 重い依存（`convx`・`kiddo`・`faer`）は `scattered` feature に閉じ込め、1D / curve / grid は軽量にビルドできる。
5. **`std` 必須**: v0.1.0 は `std` 前提。`no_std` は対象外（ロードマップ）。

---

## 2. 外部依存と技術スタック (Technology Stack)

| 役割                 | 採用クレート     | 用途・選定理由                                                                                                 |
| :------------------- | :--------------- | :------------------------------------------------------------------------------------------------------------- |
| **幾何計算エンジン** | **`convx`**      | Pure Rust の $N$ 次元 Delaunay 単体分割・隣接トポロジー。`scattered` のみ。                                    |
| **線形代数エンジン** | **`faer`**       | RBF の密行列 LU 分解のみ。`default-features = false, features = ["std"]`、版は `convx` と揃える（現在 0.24）。 |
| **空間近傍探索**     | **`kiddo`**      | k-d tree（const generic 次元）。`NearestND`・Jump・RBF 局所近傍。`scattered` のみ。                            |
| **スカラー抽象化**   | **`num-traits`** | 1D / curve / grid の `Float` ジェネリクス（f32 / f64）。                                                       |
| **並列処理 (opt)**   | **`rayon`**      | 並列バッチ評価。                                                                                               |
| **配列連携 (opt)**   | **`ndarray`**    | データ受け渡し。                                                                                               |
| **直列化 (opt)**     | **`serde`**      | 補間器の保存・復元（§6.5）。                                                                                   |
| **エラー定義**       | **`thiserror`**  | エラー型。                                                                                                     |

* 三重対角系（CubicSpline 等）は **O(n) の Thomas 法を自前実装**する（周期境界は Sherman–Morrison による巡回三重対角）。`faer` は使わない。
* `convx` は現状 crates.io 未公開（ローカル 0.0.0、ライセンス未設定）。`interpx` の公開前に `MIT OR Apache-2.0` を設定して公開する。開発中は `version` 併記の path 依存を使う。
* MSRV は `convx`・`kiddo` に合わせ **1.89**。

---

## 3. 型システムとインターフェース (Type System & APIs)

### 3.1 スカラー型と次元数
* **1D / curve / grid**: `T: num_traits::Float` によるジェネリクス（`f32` / `f64`）。
* **scattered**: **`f64` 固定**（`convx` が f64 専用のため）。f32 入力は呼び出し側が変換する。暗黙キャストは行わない。
* 次元は **`const D: usize` のみ**（grid・scattered とも）。動的次元ラッパーは v0.1.0 の対象外（ロードマップ）。`kiddo` は低次元（目安 10 未満）向けであり、doc に推奨次元を明記する。

### 3.2 共通トレイト `Interpolator`
```rust
pub trait Interpolator {
    /// クエリ点の型（1D: T, ND: [T; D] または [T]）
    type Point: ?Sized;
    /// 1 出力の型（1D: T, BSplineCurve<T, D>: [T; D]）
    type Output;

    /// 出力列数（スカラー出力は 1）。SciPy の values (N, m) の m に対応。
    fn n_outputs(&self) -> usize { 1 }

    /// 領域外（Extrapolate::None・凸包外）または NaN クエリで None。
    fn eval(&self, point: &Self::Point) -> Option<Self::Output>;

    fn eval_or(&self, point: &Self::Point, fill: Self::Output) -> Self::Output;

    fn eval_batch(&self, points: &[&Self::Point]) -> Vec<Option<Self::Output>>;

    /// 領域外は `fill` で埋める。長さ不一致は Err(LengthMismatch)。panic しない。
    fn eval_batch_into(
        &self,
        points: &[Self::Point],
        out: &mut [Self::Output],
        fill: Self::Output,
    ) -> Result<(), InterpolateError>
    where
        Self::Point: Sized;

    /// 多出力（動的 m）用。`out.len() != n_outputs()` なら Err。
    /// 領域外は Ok(false) を返し `out` には触れない。
    fn eval_into(&self, point: &Self::Point, out: &mut [Self::Scalar]) -> Result<bool, InterpolateError>;
}
```
* 上記は仕様の骨子であり、`Scalar` 関連型・`eval_batch` の引数形などの詳細は M1 の実装時に確定する（意味論は変えない）。
* SciPy のコンストラクタ引数 `fill_value` は持たない。呼び出し時の `fill` に統一する（NaN で埋めたい場合は `T::nan()` を渡す）。
* const 次元の高速パスは各構造体の固有メソッドで提供する。

**並列版**（`rayon` feature）: `eval_batch_par` / `eval_batch_into_par` を別メソッドとして提供する。

### 3.3 外挿・境界モード
```rust
pub enum Extrapolate { None, Clamp, Extend, Periodic }   // 1D / curve / grid
pub enum BoundaryCondition {                              // CubicSpline / BSpline1D
    NotAKnot,   // 既定（SciPy に合わせる）
    Natural, Clamped, Periodic,
    Derivative { /* 端点・階数・値 */ },
}
```
* `Extend` は端の多項式をそのまま延長する（SciPy の `extrapolate=True` 相当）。
* 散乱系は外挿モードを持たず、領域外（凸包外）は `None`（`fill` で置換）のみ。

### 3.4 所有権・スレッド安全性
* すべての補間器は **owned** で、ライフタイム引数を持たない。`T: 'static` の制約は課さない。
* コンストラクタは `::new(&[...])`（借用から構築）と `::from_owned(...)`（ゼロコピー移動）を提供する。
* **全補間器が `Send + Sync`** であることを公開契約とする。`Cell` / `RefCell` / `Mutex` 等の共有可変状態は持たず、`&self` での評価は純粋とする。
* Hinted search（直前区間を覚える最適化）の状態は、呼び出し側が持つ `Hint` 値かバッチ内ローカルに置く。

### 3.5 入力検証
* 入力データ（座標・値）に NaN / Inf があれば構築時に `Err(NonFinite)`。
* クエリ点が NaN の場合は panic せず `eval` は `None`、`eval_or` は `fill` を返す。

---

## 4. エラーハンドリング仕様 (Error Handling)

```rust
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum InterpolateError {
    #[error("insufficient points: required {required}, got {actual}")]
    InsufficientPoints { required: usize, actual: usize },

    #[error("dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: usize, actual: usize },   // 動的次元の場面のみ

    #[error("length mismatch: {left} vs {right}")]
    LengthMismatch { left: usize, right: usize },

    #[error("coordinates are not strictly increasing")]
    NotSorted,

    #[error("duplicate points at indices {first} and {second}")]
    DuplicatePoints { first: usize, second: usize },

    #[error("non-finite value at index {index}")]
    NonFinite { index: usize },

    #[error("degenerate geometry: {0:?}")]
    DegenerateGeometry(DegeneracyKind),

    #[error("invalid parameter `{name}`: {reason:?}")]
    InvalidParameter { name: &'static str, reason: ParamReason },

    #[error("solver failed: {0:?}")]
    SolverFailed(SolverError),
}

#[non_exhaustive] pub enum DegeneracyKind { LowerDimensional { affine_dim: usize }, NumericalLimit }
#[non_exhaustive] pub enum SolverError { Singular, /* ... */ }
#[non_exhaustive] pub enum ParamReason { /* 範囲外・必須欠落・最小次数未満 等 */ }
```
* **`convx::ConvexHullError` は公開しない**（`convx` の版と SemVer を結合させないため）。内部で次のように変換する（`From` は `pub(crate)`）:
  * `InsufficientPoints` → `InsufficientPoints`
  * `NonFiniteCoordinate{index}` → `NonFinite{index}`
  * `LengthMismatch` → `LengthMismatch`
  * `DegenerateDimension{actual_dim,..}` → `DegenerateGeometry(LowerDimensional{affine_dim})`
  * `ExactEvaluationExhausted` → `DegenerateGeometry(NumericalLimit)`
  * `TooManyPoints` → `InvalidParameter`
* 共球点群はエラーではない（`convx` が単体分割して返す）。縮退エラーは「全点が低次元部分空間上」「統合後の相異なる点が D+1 未満」のみ。

---

## 5. モジュール構造と SciPy カバレッジ (Coverage Matrix)

### 5.1 モジュール構造
```text
interpx/
├── oned/             # 1次元補間
│   ├── linear        # numpy.interp
│   ├── cubic_spline  # CubicSpline
│   ├── pchip         # PchipInterpolator
│   ├── akima         # Akima1DInterpolator (makima 含む)
│   ├── hermite       # CubicHermiteSpline
│   ├── bspline       # make_interp_spline (任意次数 k)
│   └── nearest       # Nearest1D (nearest / previous / next)
├── curve/            # N次元パラメトリック曲線
│   └── nd_curve      # BSplineCurve
├── grid/             # N次元直交グリッド補間
│   └── regular_grid  # RegularGridInterpolator
└── scattered/        # N次元不規則散乱データ補間
    ├── nearest_nd    # NearestNDInterpolator (kiddo)
    ├── linear_nd     # LinearNDInterpolator (convx + Jump-and-Walk)
    ├── clough_tocher # CloughTocher2DInterpolator
    └── rbf           # RBFInterpolator (faer + kiddo)
```
（モジュール名は Rust の識別子制約を避け `oned` とする。公開パスの最終名は M1 で確定。）

### 5.2 網羅対応表

| 分類              | SciPy 1.18                                      | `interpx` 実装構造体                          | 内部コアエンジン                       |
| :---------------- | :---------------------------------------------- | :-------------------------------------------- | :------------------------------------- |
| **1D**            | `numpy.interp`                                  | `Linear1D`                                    | `PPoly`                                |
| **1D**            | `CubicSpline`                                   | `CubicSpline`                                 | Thomas 法 + `PPoly`                    |
| **1D**            | `PchipInterpolator`                             | `PchipInterpolator`                           | `PPoly`（単調性制限導関数）            |
| **1D**            | `Akima1DInterpolator`（`method='makima'` 含む） | `Akima1D`                                     | `PPoly`                                |
| **1D**            | `CubicHermiteSpline`                            | `CubicHermiteSpline`                          | `PPoly`                                |
| **1D**            | `make_interp_spline`（任意次数 $k$）            | `BSpline1D`                                   | `BSpline`（de Boor）                   |
| **1D**            | `interp1d` の nearest / previous / next         | `Nearest1D`                                   | 二分探索                               |
| **1D**            | `PPoly` / `BSpline` の導関数・積分              | `derivative` / `antiderivative` / `integrate` | `PPoly` ⇄ `BSpline` 相互変換           |
| **N-D Curve**     | `make_interp_spline`（N-D y）                   | `BSplineCurve<T, D>`                          | `BSpline`（ベクトル制御点）            |
| **N-D Grid**      | `RegularGridInterpolator` `nearest`             | `RegularGridInterpolator`                     | インデックス四捨五入                   |
| **N-D Grid**      | `linear`                                        | 同上                                          | テンソル積線形                         |
| **N-D Grid**      | `slinear` / `cubic` / `quintic`                 | 同上                                          | 構築時テンソル積 B-spline 係数（§6.4） |
| **N-D Grid**      | `pchip`                                         | 同上                                          | 軸ごと PCHIP 遅延畳み込み（§6.4）      |
| **N-D Scattered** | `NearestNDInterpolator`                         | `NearestNDInterpolator<D>`                    | `kiddo`                                |
| **N-D Scattered** | `LinearNDInterpolator`                          | `LinearNDInterpolator<D>`                     | `convx` + Jump-and-Walk                |
| **N-D Scattered** | `CloughTocher2DInterpolator`                    | `CloughTocher2DInterpolator`                  | `convx` + 勾配推定 + Bernstein-Bézier  |
| **N-D Scattered** | `RBFInterpolator`                               | `RBFInterpolator<D>`                          | `faer` LU + `kiddo`                    |

* `interp1d` は SciPy では旧API扱い。`kind='quadratic'` / `'slinear'` は `BSpline1D` で代替する。`interp2d` は SciPy で削除済みのため対象外。

### 5.3 対象外（v0.1.0 では実装しない。§9 ロードマップ）
* FITPACK 系: `UnivariateSpline`、`splrep`、`splprep`、`RectBivariateSpline` ほか FITPACK 2D 系
* `BarycentricInterpolator`、`KroghInterpolator`
* `NdPPoly`
* `roots` / `solve`（`PPoly` / `BSpline`）
* グリッド・散乱での導関数
* 動的次元ラッパー、`no_std`、RBF の疑似逆行列フォールバック

---

## 6. コアエンジン内部アルゴリズム詳細 (Internal Algorithms)

### 6.1 スプライン基盤（`PPoly` と `BSpline`）
* **`PPoly`**: 係数を **区間ごとに連続したフラット配列**（`[区間][次数][出力列]`、出力列が最内）で保持する。SciPy の `(k+1)×(N-1)` 配置とは異なるが、評価の局所性とベクトル化のためである。`eval` は二分探索（または Hint）で区間を特定し、ホーナー法で評価する。
* **`BSpline`**: ノットベクトルと制御点を保持し、de Boor で評価する。
* `PPoly` ⇄ `BSpline` の相互変換（`from_bspline` / `to_ppoly`）、`derivative(nu)` / `antiderivative` / `integrate(a, b)` / `eval_derivative(x, nu)` を提供する。
* 三重対角系は Thomas 法、周期境界は Sherman–Morrison。

### 6.2 散乱データ幾何探索（`convx` 連携）
* **Delaunay 分割**: `[[f64; D]]` をコピーなしで `&[f64]`（row-major）にフラット化して `convx` に渡す。
* **重複点**: `convx` の `representative()` は同一座標を最小 index に統合する。`NearestND` / `LinearND` / `CloughTocher2D` は**エラーにせず受理**し、値は最小 index の点のものを使う（doc に明記）。`RBFInterpolator` は `Err(DuplicatePoints)`。
* **Point Location (Jump-and-Walk)**:
  1. **Jump**: `kiddo` で最近傍入力頂点を求め、その頂点を含む単体を開始点とする（頂点→単体の逆引きを `interpx` が構築・保持）。
  2. **Walk**: `neighbors` を辿ってクエリ点方向へ移動。ステップ上限は `n_simplices`。上限到達・巡回検知時は全単体の総当たりにフォールバックする。
  3. **重心座標**: 内部判定は `w_i >= -eps`、`eps = 100 * f64::EPSILON`。総当たりでは `eps_broad = sqrt(f64::EPSILON)`。境界上の点は許容誤差内なら内部扱い（SciPy に合わせる）。
  4. 総当たりでも見つからなければ `None`（凸包外）。
* **境界センチネル**: `NO_NEIGHBOR = u32::MAX`。`convx` は公開定数を持たないため、`interpx` 内の一箇所に定義し、`convx` 側での公開定数化を依頼する。
* `convx` の `simplices` 順序はバージョン間で保証されないため、再三角形分割に依存しない（§6.5）。

### 6.3 `CloughTocher2DInterpolator` と `RBFInterpolator`
**CloughTocher2D**
* SciPy の `estimate_gradients_2d_global` を再現する反復勾配推定（`tol = 1e-6`、`maxiter = 400`、いずれもパラメータ化）。
* 未収束でも `Err` にせず結果を使う。反復回数・`converged` を取得メソッドで返す。
* 数値一致テストの許容誤差は反復法のため緩める（§8.2）。

**RBF**
* 連立一次方程式:
  $$
  \begin{bmatrix} A & P \\ P^T & 0 \end{bmatrix}
  \begin{bmatrix} w \\ c \end{bmatrix}
  =
  \begin{bmatrix} y \\ 0 \end{bmatrix},
  \quad A_{ij} = \phi(\|x_i - x_j\|) + \lambda \delta_{ij}
  $$
  （$\lambda$ は `smoothing`。符号は SciPy 1.18 の実装に合わせて M3 で確定）
* **カーネル 8 種**（`RbfKernel` enum）: `linear`、`thin_plate_spline`、`cubic`、`quintic`、`gaussian`、`multiquadric`、`inverse_multiquadric`、`inverse_quadratic`。
* `epsilon` は gaussian・multiquadric・inverse 系で必須。必要なのに未指定なら `Err(InvalidParameter)`。
* `degree` の既定はカーネルごとの最小次数。最小次数未満の指定は `Err(InvalidParameter)`（SciPy は警告のみ）。多項式項は座標をシフト・スケーリングしてから構成する。
* 既定値: `kernel = thin_plate_spline`、`smoothing = 0`、`neighbors = None`。
* **ソルバー**: 係数行列は鞍点系で正定値ではないため、コレスキーは使わず `faer` の `partial_piv_lu` を用いる。特異なら `Err(SolverFailed(Singular))`。疑似逆行列フォールバックはロードマップ。
* **メモリ**: 全域モードは O(N²)（N = 1万点で約 800 MB）。doc に明記し、大きい N では `neighbors` 指定を推奨する。上限は設けないが、確保サイズの `usize` オーバーフローは検査して `Err`。
* `neighbors` 指定時はクエリごとに `kiddo` で近傍 $k$ 点を抽出し、局所問題を解く。

### 6.4 `RegularGridInterpolator`
* `values` は **C 順（row-major、最後の軸が最速）** の平坦な `Vec<T>` と `shape: [usize; D]` で受け取る（`ndarray` 標準レイアウト・numpy と一致）。
* `method`: `nearest` / `linear` / `slinear` / `cubic` / `quintic` / `pchip`。
* **`slinear` / `cubic` / `quintic`**: 軸ごとの補間演算子は可換なので、構築時に全軸のテンソル積 B-spline 係数を一度だけ解く。評価は局所 $(k+1)^D$ 点のみ。メモリはグリッド点数に比例。
* **`pchip`**: 値に対して非線形で軸順が結果に影響するため、SciPy と同じ軸順で、評価時に局所ステンシルのみを使って遅延畳み込みする。
* 外挿は `Extrapolate`（§3.3）。

### 6.5 直列化（`serde` feature）
* 頂点配列・隣接配列・値・パラメータをすべて直列化する（再三角形分割はしない）。`convx` のバージョンと保存データを結合させないため。
* 保存形式にフォーマット版数を持たせる。読み込み時に長さ・index 範囲を検証し、不正なら `Err`。保存済みトポロジーの Delaunay 性までは再検証しない。

---

## 7. Cargo Features 構成 (Crate Features)

```toml
[features]
default = ["scattered"]

# 1D / curve / grid は常時ビルド（重い依存なし）
scattered = ["dep:convx", "dep:kiddo", "dep:faer"]   # 外せば散乱系と重い依存が消える

# opt-in
rayon   = ["dep:rayon"]
ndarray = ["dep:ndarray"]
serde   = ["dep:serde"]
```
* `std` feature は設けない（常に `std`）。`convx` / `faer` / `kiddo` を個別の公開 feature としては出さない。
* `faer` は `default-features = false, features = ["std"]`。

---

## 8. 品質保証 (Testing, Precision, Benchmarks)

### 8.1 テスト
* **golden データ**: SciPy 1.18 を版固定した Python スクリプト（`uv`）で生成し、`tests/golden/` にコミットする。CI では SciPy を実行しない。
* **`proptest` 性質テスト**: データ点での一致、PCHIP の単調性、`LinearND` による線形関数の厳密再現、CubicSpline の導関数連続性、など。

### 8.2 数値精度基準
| 対象                             | `rtol` / `atol`   |
| :------------------------------- | :---------------- |
| 直接法（f64）                    | `1e-12` / `1e-12` |
| 反復法（CloughTocher 勾配）・RBF | `1e-6` 程度       |
| f32                              | `1e-5`            |

### 8.3 ベンチマーク
* workspace の別 crate `benches/compare`（`publish = false`）に置き、`interpx` 本体の依存に比較対象を混ぜない。
* **速度**: `criterion` で、SciPy 1.18（Python 計測スクリプト）と既存の Pure Rust クレートを同一データ・同一クエリで比較する。
* **メモリ**: カウンティングアロケーターでピーク確保量を測定（構築時と評価時を分離。RBF の O(N²) も実測）。
* 比較対象クレートは実装着手時に、メンテされている Pure Rust のクレートを調査して選定する。比較不能なカテゴリは「該当なし」と表に明記する。

---

## 9. 運用・マイルストーン・ロードマップ

### 9.1 運用
* **MSRV**: 1.89 / **edition**: 2021 / **ライセンス**: `MIT OR Apache-2.0`
* **CI（GitHub Actions）**: `fmt`、`clippy -D warnings`、Linux / Windows / macOS でのテスト、MSRV ジョブ、feature マトリクス（`--no-default-features` / `--all-features`）。

### 9.2 リリース方針
**段階公開は行わない。** 主要機能がすべて揃った時点で `v0.1.0` を一括公開する（公開による下位互換制約を早期に発生させないため）。以下は内部マイルストーンであり、公開バージョンではない。`convx` が公開されるまでは、`convx` 依存のない M1〜M3 を先に進める。

- [ ] **M1: Core + 1D + Curve**
  - [ ] `Interpolator` トレイト・`InterpolateError`・`Extrapolate` / `BoundaryCondition`
  - [ ] `PPoly` / `BSpline`（導関数・積分・相互変換、Thomas 法）
  - [ ] `Linear1D`、`Nearest1D`、`CubicSpline`（NotAKnot / Natural / Clamped / Periodic / Derivative）
  - [ ] `PchipInterpolator`、`Akima1D`（makima）、`CubicHermiteSpline`、`BSpline1D`
  - [ ] `BSplineCurve`
- [ ] **M2: Grid**
  - [ ] `RegularGridInterpolator`（`nearest` / `linear` / `slinear` / `cubic` / `quintic` / `pchip`）
- [ ] **M3: `convx` 非依存の散乱系**
  - [ ] `NearestNDInterpolator`（`kiddo`）
  - [ ] `RBFInterpolator`（`faer` LU、8 カーネル、多項式項、局所近傍）
- [ ] **M4: `convx` 依存の散乱系**（`convx` 公開後に着手）
  - [ ] `LinearNDInterpolator`（Delaunay + Jump-and-Walk、§6.2）
  - [ ] `CloughTocher2DInterpolator`（勾配推定 + Bernstein-Bézier）
- [ ] **M5: 統合と品質**
  - [ ] `rayon` 並列バッチ評価、`ndarray` 相互運用、`serde`
  - [ ] SciPy 1.18 golden テストスイート、`proptest`
  - [ ] 比較ベンチ（SciPy・既存 Rust クレート、速度・メモリ）
  - [ ] docs（推奨次元、RBF メモリ、重複点の扱い 等）、CI
  - [ ] **完了で `v0.1.0` を公開**

### 9.3 ロードマップ（v0.1.0 以降）
* `roots` / `solve`（`PPoly` / `BSpline`）、`NdPPoly`
* グリッド・散乱での導関数
* 動的次元ラッパー（grid・scattered）
* RBF の疑似逆行列フォールバック
* `BarycentricInterpolator`、`KroghInterpolator`
* FITPACK 系（`UnivariateSpline`、`splrep`、`splprep`、`RectBivariateSpline` ほか）
* `no_std`（`alloc`）対応の検討
