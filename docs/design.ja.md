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

### 2.1 開発用依存 (dev-dependencies)

テスト・ベンチ専用で、公開パッケージの依存には現れない。ここに名指ししたもの以外は追加しない（追加は設計変更が先。`.cursor/rules/layout.mdc`）。

| クレート         | 用途                                                                                       | 最初に必要とする行                   |
| :--------------- | :----------------------------------------------------------------------------------------- | :----------------------------------- |
| **`serde_json`** | golden JSON（`compare/goldens/*.json`）の読み込み。`tests/common` のローダーで使う。      | M1-1                                 |
| **`serde`**      | 上記ローダーで JSON を型付き構造体に読むための `derive`（`serde` が dev で必要な場合のみ）。 | M1-1                                 |
| **`proptest`**   | 性質テスト（§8.1）。                                                                       | 最初に性質テストを書く行（M1-3 以降）。 |
| **`criterion`**  | `benches/` の時間計測（§8.3）。CI では実行しない。                                         | 最初のベンチを追加する行             |
| **`stats_alloc`** | `tests/alloc.rs` のカウンティングアロケーター（§8.1、§8.3 のメモリ計測にも使う）。         | M1-2                                 |

* ゼロアロケーション検査（`tests/alloc.rs`、§8.1）のカウンティングアロケーターは、**外部クレート `stats_alloc` を使う**。`GlobalAlloc` の実装には `unsafe impl` が必要で、`unsafe_code = "forbid"` は `interpx` 自身のコード（テストを含む）にかかるため、自前実装はできない。`unsafe` は依存側に置く（`gprx` も同じクレートを使う）。
* 各依存のバージョンと features は、依存を追加する実装行で決める。

---

## 3. 型システムとインターフェース (Type System & APIs)

### 3.1 スカラー型・値の型・次元数
* **スカラー `Scalar`**: 封印したトレイトで、`f32` と `f64` にだけ実装する。`num_traits::Float + Debug + Send + Sync` を含む。補間器の型引数 `T` の境界はこれ 1 つにする。
* **1D / curve / grid**: `T: Scalar` によるジェネリクス（`f32` / `f64`）。
* **scattered**: **`f64` 固定**（`convx` が f64 専用のため）。f32 入力は呼び出し側が変換する。暗黙キャストは行わない。
* **値の型 `Value`**: 補間する値（SciPy の `values` の 1 行）の型。封印したトレイトで、スカラー `T` と固定長配列 `[T; M]` に実装する。必要な演算はゼロ、加算つきスカラー倍、有限性の判定だけである。

  ```rust
  pub trait Value<T: Scalar>: Copy + Debug + Send + Sync {
      fn zero() -> Self;
      /// self + k * x
      fn axpy(self, k: T, x: Self) -> Self;
      fn all_finite(&self) -> bool;
  }
  ```

  * 多出力（SciPy の `values` が `(N, m)`）は、値の型を `[T; M]` にして表す。例: `LinearNDInterpolator<const D: usize, V = f64>`、`RegularGridInterpolator<T, const D: usize, V = T>`、`BSplineCurve<T, D>`（値の型は `[T; D]`）。
  * `[T; M]` は連続領域なので、`PPoly` の係数レイアウト `[区間][次数][出力列]`（§6.1）と一致する。
  * 出力列数 `m` を実行時に決める形は v0.1.0 の対象外（ロードマップ A-4。動的次元と同じ扱い）。
* 次元は **`const D: usize` のみ**（grid・scattered とも）。動的次元ラッパーは v0.1.0 の対象外（ロードマップ）。`kiddo` は低次元（目安 10 未満）向けであり、doc に推奨次元を明記する。

### 3.2 共通トレイト `Interpolator`
```rust
pub trait Interpolator {
    /// クエリ点の型。1D は T、N-D は [T; D]。
    type Point;
    /// 出力の型。T、または [T; M]。
    type Output;

    /// 1 点を評価する。領域外（`Extrapolate::Forbid`・凸包外）と NaN クエリは None。
    fn eval(&self, point: &Self::Point) -> Option<Self::Output>;

    /// 領域外と NaN クエリで `fill` を返す。
    fn eval_or(&self, point: &Self::Point, fill: Self::Output) -> Self::Output {
        self.eval(point).unwrap_or(fill)
    }

    /// 結果を新しい Vec に入れて返す（確保する）。
    fn eval_batch(&self, points: &[Self::Point]) -> Vec<Option<Self::Output>> {
        points.iter().map(|p| self.eval(p)).collect()
    }

    /// 結果を `out` に書く（確保しない）。領域外と NaN クエリは `fill` で埋める。
    /// `points.len() != out.len()` は Err(LengthMismatch) で、`out` には触れない。
    fn eval_batch_into(
        &self,
        points: &[Self::Point],
        out: &mut [Self::Output],
        fill: Self::Output,
    ) -> Result<(), InterpolateError>
    where
        Self::Output: Clone;
}
```
* 必須メソッドは `eval` だけで、残りは既定実装を持つ。各補間器は、速くできる場合に `eval_batch_into` を上書きする（区間探索の `Hint` をバッチ内で使い回す、など）。
* `Point` は `Sized`（v0.1.0 は `const D` のみ）。トレイトは封印しない（利用者が自分の補間器に実装できる）。
* SciPy のコンストラクタ引数 `fill_value` は持たない。呼び出し時の `fill` に統一する（NaN で埋めたい場合は `T::nan()` を渡す）。
* `eval`・`eval_or`・`eval_batch_into` は確保しない（§8.1 のゼロアロケーション検査）。`eval_batch` は `Vec` を返すので確保する。

**並列版**（`rayon` feature）: `eval_batch_par` / `eval_batch_into_par` を別メソッドとして提供する。

### 3.3 外挿・境界モード
```rust
#[non_exhaustive]
pub enum Extrapolate { Forbid, Clamp, Extend, Periodic }   // 1D / curve / grid

pub enum BoundaryCondition<T> {                            // CubicSpline / BSpline1D
    NotAKnot,                                 // 既定（SciPy に合わせる）
    Periodic,
    Ends(EndCondition<T>, EndCondition<T>),   // 左端、右端
}
pub enum EndCondition<T> { NotAKnot, FirstDerivative(T), SecondDerivative(T) }
```
* **`Extrapolate`**: `Forbid` は領域外で `None`。`Clamp` は端の値を返す（`numpy.interp` の挙動）。`Extend` は端の多項式をそのまま延長する（SciPy の `extrapolate=True`）。`Periodic` は領域に折り返す。
* **既定値は型ごとに SciPy 1.18 に合わせる。**

  | 型                                                                      | 既定        |
  | :---------------------------------------------------------------------- | :---------- |
  | `Linear1D`                                                              | `Clamp`     |
  | `CubicSpline`、`PchipInterpolator`、`CubicHermiteSpline`、`BSpline1D`、`BSplineCurve` | `Extend`    |
  | `CubicSpline`・`BSpline1D` で境界条件が `Periodic` のとき               | `Periodic`  |
  | `Akima1D`、`Nearest1D`、`RegularGridInterpolator`                       | `Forbid`    |

* **`BoundaryCondition`**: 端ごとに条件を持ち、左右で別の条件を指定できる（SciPy の `bc_type=((1, a), (2, b))`）。`Natural` は両端 `SecondDerivative(0)`、`Clamped` は両端 `FirstDerivative(0)` で、`BoundaryCondition::natural()` / `clamped()` で作る。周期は両端まとめての指定で、片側だけ周期という状態は型で表せない。
* 散乱系は外挿モードを持たず、領域外（凸包外）は `None`（`fill` で置換）のみ。

### 3.4 所有権・スレッド安全性
* すべての補間器は **owned** で、ライフタイム引数を持たない。`T: 'static` の制約は課さない。
* コンストラクタは `::new(&[...])`（借用から構築）と `::from_owned(...)`（ゼロコピー移動）を提供する。
* **全補間器が `Send + Sync`** であることを公開契約とする（`Scalar` と `Value` の境界が保証する）。`Cell` / `RefCell` / `Mutex` 等の共有可変状態は持たず、`&self` での評価は純粋とする。
* Hinted search（直前区間を覚える最適化）の状態は、呼び出し側が持つ `Hint` 値かバッチ内ローカルに置く。

### 3.5 入力検証
* 入力データ（座標・値）に NaN / Inf があれば構築時に `Err(NonFinite)`。
* クエリ点が NaN の場合は panic せず `eval` は `None`、`eval_or` は `fill` を返す。
* **1D の座標は狭義単調増加のみ**受け付ける（SciPy の `CubicSpline` と同じ）。降順の軸を受け付けるのはグリッドだけである。
* **検証の順序は固定する。** 複数の誤りがある入力は、次の順で最初に当たったエラーを返す。
  1. `LengthMismatch`（座標と値の長さ）
  2. `InsufficientPoints`
  3. `NonFinite`（座標を先に、次に値。最初の index）
  4. `DuplicatePoints`（隣り合う座標が等しい）と `NotSorted`（減少）のうち、先に現れたほう

---

## 4. エラーハンドリング仕様 (Error Handling)

```rust
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum InterpolateError {
    #[error("insufficient points: required {required}, got {actual}")]
    InsufficientPoints { required: usize, actual: usize },

    #[error("length mismatch: expected {expected}, got {actual}")]
    LengthMismatch { expected: usize, actual: usize },

    #[error("non-finite value at index {index}")]
    NonFinite { index: usize },

    #[error("coordinates are not strictly increasing at index {index}")]
    NotSorted { index: usize },

    #[error("duplicate points at indices {first} and {second}")]
    DuplicatePoints { first: usize, second: usize },

    #[error("periodic boundary needs equal values at both ends")]
    NotPeriodic,
}
```
* `#[non_exhaustive]` なので、variant は最初に必要になる行で追加する（使う前に定義しない）。上は M1（1D・curve）が使うものである。
* 後の行で追加する variant:

  | variant                                                    | 追加する行         | 用途                                                   |
  | :--------------------------------------------------------- | :----------------- | :----------------------------------------------------- |
  | `InvalidParameter { name: &'static str, reason: ParamReason }` | M1-4（`BSpline`） | 引数が範囲外。M1-4 は `ParamReason::TooLarge { max: usize }`（次数の上限）を入れる。必須欠落・最小次数未満などの `ParamReason` は、それを返す行で足す |
  | `SolverFailed(SolverError)`                                | M3（RBF）          | 特異な連立方程式（`SolverError::Singular`）            |
  | `DegenerateGeometry(DegeneracyKind)`                       | M4（`LinearND`）   | `LowerDimensional { affine_dim }`、`NumericalLimit`    |

* `NonFinite { index }` の `index` は、座標に非有限値があればその座標の位置で、座標がすべて有限のときに限り値の位置である（検証は座標が先。§3.5）。variant からはどちらの配列かを区別できない。
* 文字列のペイロードは持たない（パターンマッチできる型付きの値だけ）。`ParamReason`・`SolverError`・`DegeneracyKind` も `#[non_exhaustive]`。
* 次元の不一致を表す variant は持たない（`const D` のみなので、型で表現できない）。
* `NotSorted { index }` の `index` は、直前の座標より小さい最初の要素。`LengthMismatch` の `expected` は座標（またはクエリ）の長さ、`actual` は値（または出力バッファ）の長さ。
* **`convx::ConvexHullError` は公開しない**（`convx` の版と SemVer を結合させないため）。M4 で、内部で次のように変換する（`From` は `pub(crate)`）:
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
  * **型**: `PPoly<T: Scalar, V: Value<T> = T>`。次数 $k$ は実行時の値で持つ。区切り点 $x_0 < x_1 < \dots < x_n$ と、区間ごとに $k+1$ 個の係数を持つ。
  * **係数の並び**: 区間 $i$ の値は $\sum_{j=0}^{k} c_{i,j}\,(x - x_i)^{k-j}$ で、区間内は**高次の項が先**。ホーナー法がメモリを前から順に読める。SciPy の `c[power][interval]` も高次が先なので、同じ数を区間ごとに並べ替えただけである。
  * **構築**: `PPoly::new(&[T], &[V], degree)`（複製）と `from_owned(Vec<T>, Vec<V>, degree)`（移動）。次数は引数で明示し、係数の長さから逆算しない。検証は §3.5 の順で、係数の個数が `区間数 × (degree + 1)` と違えば `LengthMismatch`、区切り点が 2 個未満なら `InsufficientPoints`。
  * **区間の帰属**: 区間は左閉右開 $[x_i, x_{i+1})$ で、最後の区間だけ右端を含む（SciPy と同じ）。
  * **領域外**: `Extrapolate` の 4 モード。既定は `Extend`（SciPy の `PPoly` と同じ）で、`with_extrapolate` で変える。`Periodic` は $x_0 + ((x - x_0) \bmod (x_n - x_0))$ に折り返す。NaN のクエリは `None`。無限大のクエリは `Forbid` と `Periodic` で `None`、`Clamp` は端の値、`Extend` は端の多項式をそのまま評価した値（一般に無限大か NaN）。
  * **`Hint`**: `eval_hinted(x, &mut Hint)` は、ヒントの区間とその両隣を先に調べ、外れたら二分探索に戻る。結果は `eval` とビット同一で、変わるのは時間だけである。
  * **ヒントが効く条件**: ヒントつきの探索は直前の点の区間から始めるので、各点が直前の点の結果を待つ。近い点が続けば待ちは短く、二分探索を省ける分だけ速い（昇順 10000 点で約 2 倍）。順不同の点では待ちが二分探索まるごとになり、互いを待たない単純な探索より約 2 倍遅い（`benches/ppoly.rs` の実測）。`eval_hinted` は順に並んだクエリ用で、順不同なら `eval` を使う。
  * **バッチ**: `eval_batch_into` は 64 点ずつに区切り、区切りごとにヒントを使うかを決める。使うのは、区切りが単調（昇順か降順）で、またぐ区間の数が点数以下のときだけである。16 点未満の区切りは判定せず単純に探索する（判定の固定コストのほうが大きい）。どの経路でも各点の区間は `eval` と同じで、結果はビット同一である。
* **`BSpline`**: ノット列と制御点を保持し、非ゼロの基底関数と制御点の線形結合で評価する。
  * **型**: `BSpline<T: Scalar, V: Value<T> = T>`。次数 $k$ は実行時の値で持つ。制御点 $c_0, \dots, c_{n-1}$ と、ノット $t_0 \le t_1 \le \dots \le t_{n+k}$（重複を許す）を持つ。
  * **構築**: `BSpline::new(&[T], &[V], degree)`（複製）と `from_owned(Vec<T>, Vec<V>, degree)`（移動）。検証は次の順で、最初に当たったエラーを返す。
    1. `InvalidParameter { name: "degree", reason: TooLarge { max: 15 } }`（次数が上限を超える）
    2. `LengthMismatch`（ノットの個数が `制御点の個数 + degree + 1` と違う）
    3. `InsufficientPoints`（制御点が `degree + 1` 個未満）
    4. `NonFinite`（ノットを先に、次に制御点）
    5. `NotSorted`（ノットが減少する）
    6. `DuplicatePoints { first: k, second: n }`（定義域の幅がゼロ、すなわち $t_k = t_n$）
  * **定義域と区間の帰属**: 定義域は $[t_k, t_n]$。点 $x$ の区間は $t_\ell \le x < t_{\ell+1}$ を満たす幅が正のノット区間で、定義域の右端 $t_n$ は最後の幅が正の区間に含める（`PPoly` と同じ左閉右開、最後だけ閉）。重複ノットが作る幅ゼロの区間は、どの点の区間にもならない。
  * **評価**: 区間 $\ell$ で非ゼロの基底関数 $B_{\ell-k}, \dots, B_\ell$ の値（$k+1$ 個）を de Boor–Cox の漸化式でスタック上の配列に求め、$\sum_j B_{\ell-k+j}(x)\,c_{\ell-k+j}$ を `Value::axpy` で足し合わせる。評価する区間 $\ell$ は幅が正なので、漸化式が割る数（$t_{\ell+n} - t_{\ell+n-j}$）はどれも区間 $\ell$ をまたぎ、ノットが重複していても正である。SciPy の実装（`_deBoor_D`）にある「等しいノットの項を飛ばす」分岐は要らない。
  * **漸化式の形**: 値は、$x$ から両側のノットまでの距離を保持して値をその場で更新する三角形の形（Piegl–Tiller の A2.2）で求める。SciPy と同じ「次数ごとに作業配列を写す」形は、導関数（`from_bspline` のテイラー係数、および M1-5）にだけ使う。値の計算を後者で行うと、3 次で約 2 倍遅かった（`benches/bspline.rs` の実測。1000 点で 40.8 µs に対し 21.4 µs）。除数の逆数を先に取る形は差が誤差の範囲だったので採らない。
  * **次数の上限**: 評価は確保しない（§3.2）。作業用の配列をスタックに固定長で置くため、次数は 15 までとする。SciPy の `make_interp_spline` が実用上使うのは 5 次までである。
  * **領域外**: `Extrapolate` の 4 モードで、既定は `Extend`（SciPy の `BSpline` と同じ）。`Extend` は最初または最後の幅が正の区間の多項式を延長する。`Periodic`・NaN・無限大のクエリの扱いは `PPoly` と同じ。
  * **区間探索の共有**: 二分探索、`Hint` つきの探索、バッチで区切りごとにヒントを使うかの判定は、crate 内のモジュール `search` に置き、`PPoly` と `BSpline` が共用する。`Hint` 型も共通で、`BSpline` も `eval_hinted` を持つ。
* **`PPoly` と `BSpline` の変換**: どちらも厳密な変換（丸め誤差を除く）で、外挿モードを引き継ぐ。
  * **`PPoly::from_bspline(&BSpline) -> Result<PPoly, _>`**（SciPy の `PPoly.from_spline`）: 区切り点は定義域 $[t_k, t_n]$ 内の相異なるノット。区間 $[u_j, u_{j+1}]$ の係数は、スプラインの $m$ 階導関数の $u_j$ での値を $m!$ で割ったもの（$u_j$ まわりのテイラー係数）。SciPy は全ノットを区切り点に残すので重複ノットが幅ゼロの区間になるが、`PPoly` の区切り点は狭義単調増加なので、幅ゼロの区間は作らない。
  * **`BSpline::from_ppoly(&PPoly) -> Result<BSpline, _>`**: 各区切り点に多重度 $k+1$ のノットを置き、区間ごとにべき基底の係数を Bernstein 基底の係数（Bézier の制御点）に直す。区間どうしの連続性を仮定しないので、不連続な区分多項式でも厳密に表せる。制御点は `区間数 × (k + 1)` 個になる。連続性を調べてノットの多重度を減らす処理は v0.1.0 の対象外（ロードマップ A-9）。次数が 15 を超える `PPoly` は `InvalidParameter` を返す。
  * 往復（`PPoly` → `BSpline` → `PPoly`、およびその逆）で、評価値が許容誤差内で一致することをテストする。
* `derivative(nu)` / `antiderivative` / `integrate(a, b)` / `eval_derivative(x, nu)` を、`PPoly` と `BSpline` に提供する（M1-5）。
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
* **SciPy 1.18 の既定値との差（`cubic` / `quintic`）**: SciPy はこの 2 メソッドの係数を、既定では反復解法（`gcrotmk`）で求める。その結果はテンソル積スプラインの厳密解から、等間隔の 2D グリッドで約 $10^{-5}$、非等間隔や 3D のグリッドでは最大約 $10^{-2}$ ずれる（golden 作成時の実測）。`interpx` は係数を直接法で厳密に解くため、SciPy の既定の出力とはこの分だけ一致しない。一致する相手は、SciPy で `solver=scipy.sparse.linalg.spsolve` を指定した結果、および `cubic_legacy` / `quintic_legacy` である（両者は $3 \times 10^{-13}$ 以内で一致する）。golden の `expected` は直接解法の値とし、既定の値は `expected_scipy_default` として併記する（§8.4.1）。この差は rustdoc にも書く。
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
* **golden データ**: `compare/` の Python スクリプト（`compare/generate_*.py`、`uv` で SciPy 1.18 を `pyproject.toml` / `uv.lock` に固定）で生成し、`compare/goldens/*.json` にコミットする。`cargo test` は読むだけで Python を呼ばず、CI でも SciPy を実行しない。再生成は `just gen-goldens`（M1-1 で追加）。`compare/goldens/` は公開パッケージに含める（`Cargo.toml` の `include`）。
* **ゼロアロケーション検査**: `tests/alloc.rs` がカウンティングアロケーターで、`eval` / `eval_or` / `eval_batch_into` の確保回数が 0 であることを検査する（`just test` に含む）。並列版 `eval_batch_par` は逐次版とビット同一であることをテストする。
* **`proptest` 性質テスト**: データ点での一致、PCHIP の単調性、`LinearND` による線形関数の厳密再現、CubicSpline の導関数連続性、など。

### 8.2 数値精度基準
| 対象                             | `rtol` / `atol`   |
| :------------------------------- | :---------------- |
| 直接法（f64）                    | `1e-12` / `1e-12` |
| 反復法（CloughTocher 勾配）・RBF | `1e-6` 程度       |
| f32                              | `1e-5`            |

### 8.3 ベンチマーク
* ルート crate は単一 crate（workspace なし）。比較用の資材はルート直下の `compare/`（`pyproject.toml`、`uv.lock`、`generate_*.py`、`goldens/`、`perf/`）に置く。`compare/perf/` の Rust ランナーは、他クレートに依存する独立した Cargo パッケージ（自前の `[workspace]`）とし、`interpx` 本体の依存に比較対象を混ぜない。
* 比較実行は手動で、CI のゲートではない（結果は `compare/perf/out/`、gitignore 対象）。`criterion` ベンチは `benches/` に置く。
* **速度**: `criterion` で、SciPy 1.18（Python 計測スクリプト）と既存の Pure Rust クレートを同一データ・同一クエリで比較する。
* **メモリ**: カウンティングアロケーターでピーク確保量を測定（構築時と評価時を分離。RBF の O(N²) も実測）。
* 比較対象クレートは実装着手時に、メンテされている Pure Rust のクレートを調査して選定する。比較不能なカテゴリは「該当なし」と表に明記する。

### 8.4 検証・ベンチ用データセット

テスト（§8.1）とベンチ（§8.3）が使うデータを、目的で 4 層に分ける。どの行（`docs/roadmap.md`）も、ここに無いデータを黙って追加しない。追加は設計変更が先。

| 層               | 目的                             | 基準・規模                                         | 実行                                 |
| :--------------- | :------------------------------- | :------------------------------------------------- | :----------------------------------- |
| **A 正しさ**     | SciPy 1.18 との数値一致          | 小規模・決定的。`compare/goldens/` にコミット      | `cargo test`                         |
| **B 収束・性質** | 解析解との差、収束次数、不変条件 | 解析的に値が分かる関数。`proptest`                 | `cargo test`                         |
| **C ベンチ**     | 速度・メモリ                     | SciPy のベンチと同じパラメータ＋ Rust 単独の大規模 | 手動（`just bench`、`compare/perf`） |
| **D 敵対的**     | 退化・近接・境界での挙動         | 小規模・決定的                                     | `cargo test`                         |

#### 8.4.1 層 A: 正しさ

**Franke の 2D 関数（定義域 $[0,1]^2$）。** Franke (1979) の散乱データ補間の比較研究で使われた 6 曲面で、元は McLain (1974) の関数である。$r^2 = (x - \tfrac12)^2 + (y - \tfrac12)^2$ とおく。

| ID | 名称                 | 定義                                                                                                                                                                         |
| :- | :------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| F1 | 総合                 | $0.75\,e^{-\frac{(9x-2)^2 + (9y-2)^2}{4}} + 0.75\,e^{-\frac{(9x+1)^2}{49} - \frac{(9y+1)^2}{10}} + 0.5\,e^{-\frac{(9x-7)^2 + (9y-3)^2}{4}} - 0.2\,e^{-(9x-4)^2 - (9y-7)^2}$ |
| F2 | 崖（Cliff）          | $\dfrac{\tanh\bigl(9(y - x)\bigr) + 1}{9}$                                                                                                                                   |
| F3 | 鞍点（Saddle）       | $\dfrac{1.25 + \cos(5.4\,y)}{6 + 6\,(3x - 1)^2}$                                                                                                                             |
| F4 | 緩やかな丘（Gentle） | $\tfrac13\,e^{-\frac{81}{16} r^2}$                                                                                                                                           |
| F5 | 急な丘（Steep）      | $\tfrac13\,e^{-\frac{81}{4} r^2}$                                                                                                                                            |
| F6 | 球（Sphere）         | $\dfrac{\sqrt{64 - 81\,r^2}}{9} - \tfrac12$（半径 $\tfrac89$、中心 $(\tfrac12, \tfrac12, -\tfrac12)$ の球面。単位正方形では根号の中は常に正）                                |

* **照合の記録**: F1〜F5 は 2 つの独立した実装（UQTestFuns、および Burkardt の `test_interp_2d`（MIT））で式が一致した（F4 の係数 $\tfrac{81}{16}$ = 5.0625、F5 の $\tfrac{81}{4}$ = 20.25、いずれも因子 $\tfrac13$ を持つ）。F6 は UQTestFuns が $\tfrac19$ を欠くが、半径 $\tfrac89$ の球という記述と整合するのは $\tfrac19$ を持つ式であり、Burkardt の実装もこれと一致する。F1 の第 2 項は Burkardt の抜粋で二乗が見えないため、UQTestFuns の式を採る。Franke (1979) 原典は取得できていないため、golden を作る最初の行で原典または Renka–Brown の公開コードと再照合し、結果をこの節に追記する。
* **散乱点の配置**: Franke の 100 点・33 点と、Lawson の 25 点の節点配置を使う（25 点は Franke ではなく Lawson の配置である）。座標は Burkardt の `test_interp_2d`（MIT）から機械的に転記する。Franke の 100 点には単位正方形の外に出る点があるため、単位正方形上のクエリの一部は凸包の外になる。
* **golden の範囲**: 層 A の golden は、補間族ごとに `compare/goldens/` の `oned` / `curve` / `grid` / `scattered` / `rbf` に置く。入力のみの `datasets`（古典データ・節点配置・敵対的入力）と、解析関数の値の `functions` も置く。1D のケースは値に加えて、1 階・2 階導関数、原始関数、定積分、区分多項式の係数（または B-spline のノットと制御点）を持つ。グリッドの `cubic` / `quintic` の基準は §6.4 のとおり直接解法の値とする。
* **サイズの上限**: golden は 1 ファイル 1 MB、ディレクトリ合計 4 MB を上限とする（公開パッケージに含めるため）。

**1D の古典データ。** 値は一次資料（原論文・原典）から転記し、出典を golden の `meta` に記録する。

| データ                   | 点数 | 見るところ                                                                 | 出典                       |
| :----------------------- | :--- | :------------------------------------------------------------------------- | :------------------------- |
| Fritsch–Carlson          | 11   | 急な方向転換とまばらな点。オーバーシュートと PCHIP の単調性                | Fritsch・Carlson (1980)    |
| チタンの物性（titanium） | 49   | 間隔が不均一なノード。スプライン近似の標準データ                           | de Boor・Rice              |
| $y = x^2$ ＋近接 1 点    | 12   | 近接した $x$ で $y$ が整合しない点に対する敏感性                           | Burkardt `test_interp` p08 |
| Runge 関数               | —    | $f(x) = \dfrac{1}{1 + 25x^2}$（$x \in [-1, 1]$、等間隔ノード）の端部の振動 | Runge                      |
| ステップ                 | —    | $x < 0.5$ で 0、それ以外で 1。単調性、`CubicSpline` のオーバーシュート     | —                          |
| チャープ                 | —    | $\sin(x^2)$（$x \in [0, 10]$）。位相と高周波の歪み                         | —                          |

#### 8.4.2 層 B: 収束と性質

* **収束次数**: 滑らかな関数で $h$ を半分にして誤差を測り、$p = \log_2 (e_h / e_{h/2})$ を求める。期待値（滑らかなデータ）は次のとおり。

  | 補間                                        | 期待される次数 |
  | :------------------------------------------ | :------------- |
  | `Linear1D`、グリッド `linear`               | 2              |
  | `CubicSpline`（NotAKnot）、グリッド `cubic` | 4              |
  | `PchipInterpolator`                         | 3              |
  | `BSpline1D`（次数 $k$）                     | $k + 1$        |
  | グリッド `quintic`                          | 6              |

  上に無い補間（Akima、Hermite、`LinearND`、`CloughTocher2D`、RBF）の次数は、その行の Issue に出典付きで書く。
* **解析解**: 導関数・不定積分・定積分は、$\sin$、多項式、指数関数など解析的に求まる関数との差を測る。
* **性質**: §8.1 の `proptest` 項目。
* 関数: 1D は上の古典関数と $\sin$、指数関数。2D は Franke F1〜F6。

#### 8.4.3 層 C: ベンチ

**SciPy 1.18 のベンチ（`benchmarks/benchmarks/interpolate.py`）と同じパラメータ**で比較する。パラメータは SciPy の版を上げるたびに再確認する。

| 対象                                         | SciPy のベンチのパラメータ                                                                                                         |
| :------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------- |
| 1D（`Interpolate1d`）                        | $n \in \{10, 50, 100, 1000, 10000\}$、メソッド linear / nearest / zero / slinear / quadratic / cubic。指数減衰データ、評価点は $4n$ |
| `PPoly` / `BSpline`                          | 区間 55 の 3 次。$n \in \{10, 100, 1000, 10000\}$ 点の評価                                                                         |
| `RegularGridInterpolator`                    | 次元 $\{2, 3, 4\}$、軸あたりの最大点数 $\{10, 40, 200\}$、評価点数 $\{10, \dots, 10000\}$、軸の向き $\pm 1$                        |
| グリッド cubic / quintic（2D）               | cubic: 軸あたり $\{10, 40, 100, 200, 400\}$、quintic: $\{10, 40\}$                                                                 |
| `griddata`（`nearest` / `linear` / `cubic`） | 2D の乱数 1000 点、出力グリッド軸あたり $\{10, 100, 1000\}$                                                                        |
| `NearestNDInterpolator`                      | $n \in \{10, 100, 1000\}$、出力グリッド軸あたり $\{50, 100, 500\}$                                                                 |
| `CloughTocher2D`                             | $n \in \{10, 50, 100\}$、$50 \times 50$ 評価                                                                                       |
| `RBFInterpolator`                            | `neighbors` が `None` または 50、$N \in \{10, 100, 1000\}$、全 8 カーネル、`epsilon = 5.0`、$[-1, 1]^2$ の乱数点                   |

* **Rust 単独の大規模**: `NearestND`・`LinearND` は $N$ を $10^5$ まで、グリッドは軸あたり数百点まで。RBF の全域モードは $N \le 10^4$（メモリは $O(N^2)$）、$N$ がそれ以上は `neighbors` 指定のみ。
* 既存の Pure Rust クレートとの比較は M5-5 で行う（§8.3）。
* 構築と評価は別々に測る。バッチ評価は `eval_batch_into` と `eval_batch_into_par` を、SciPy のベクトル化呼び出しと比べる。

#### 8.4.4 次元 $D \ge 3$ の関数

* **Genz (1984) の関数族**（$x \in [0,1]^D$、$u_i = \tfrac12$）。パラメータは $a_i = 2 + i/D$（$i = 0, \dots, D-1$）に固定する。これは実装上の選択であり、Genz の難易度の正規化ではない。式は標準的な形を採ったが原典とは未照合で、最初にこの族を補間の検証に使う行で照合して記録する。

  | 族                            | 式                                                    |
  | :---------------------------- | :---------------------------------------------------- |
  | 振動（Oscillatory）           | $\cos\bigl(2\pi u_1 + \sum_i a_i x_i\bigr)$           |
  | 積ピーク（Product peak）      | $\prod_i \bigl(a_i^{-2} + (x_i - u_i)^2\bigr)^{-1}$   |
  | コーナーピーク（Corner peak） | $\bigl(1 + \sum_i a_i x_i\bigr)^{-(D+1)}$             |
  | ガウスピーク（Gaussian peak） | $\exp\bigl(-\sum_i a_i^2 (x_i - u_i)^2\bigr)$         |
  | 連続（Continuous）            | $\exp\bigl(-\sum_i a_i \lvert x_i - u_i \rvert\bigr)$ |

* **Rosenbrock**（曲がった谷）: $\sum_{i=1}^{D-1} \bigl[100\,(x_{i+1} - x_i^2)^2 + (1 - x_i)^2\bigr]$ を 1 ケースとして使う。
* 次元は $D \in \{3, 5, 8\}$（`kiddo` の推奨は低次元、目安 10 未満）。

#### 8.4.5 散乱点のサンプリング

4 種とも、シードを固定して決定的にする。

1. **Halton 列**（準乱数）: 空間充填性のある不規則な観測点。
2. **一様乱数**: 密集と疎が自然に生じる。
3. **グリッド欠損**: 正方グリッドから 30〜70% の点を除いたもの。
4. **幾何分布**: Qhull の `rbox` に相当する分布（立方体内、球面上、立方体表面など）。合わせるのは分布の種類だけで、乱数生成器は NumPy の `Generator` を使う（`convx` の `benches/sets.txt` の生成器とは一致しない）。

#### 8.4.6 層 D: 敵対的入力

期待する挙動は、ADR 0005 と 0006 に合わせる。

| ケース                 | データ                                        | 期待する挙動                                                                               |
| :--------------------- | :-------------------------------------------- | :----------------------------------------------------------------------------------------- |
| 共円・共球             | 円周・球面上の 4 点以上（2D）、5 点以上（3D） | 成功し、結果は有限。`LinearND` はアフィン関数を厳密に再現する                              |
| 共線・共面             | 同一直線・同一平面上の点群                    | `Err(DegenerateGeometry)`                                                                  |
| 極細単体               | 1 辺が極端に短い扁平な単体                    | 重心座標は有限で、許容誤差内                                                               |
| 完全な重複点           | 同一座標の点                                  | `RBFInterpolator` は `Err(DuplicatePoints)`。他は統合して受理し、値は最小 index の点のもの |
| 極近接（重複ではない） | 距離 $10^{-12}$ 程度の点対                    | `Ok` で有限値、または `Err(SolverFailed(Singular))`。panic・無限ループは起こさない         |
| 凸包外のクエリ         | 点群の外側                                    | `None`（`eval_or` は `fill`）。無限ループにならない                                        |
| NaN クエリ             | クエリに NaN                                  | `None`（`eval_or` は `fill`）                                                              |
| 非有限の入力データ     | 座標・値に NaN / Inf                          | `Err(NonFinite)`                                                                           |

#### 8.4.7 実世界データ

* v0.1.0 に含める。ただし `cargo test` には入れず、`compare/perf` の手動比較でのみ使う。
* 取得スクリプト・版固定・チェックサムを `compare/perf/` に置き、ライセンスを確認したもののみ使う。生データはコミットしない。
* 具体的なデータセット（候補: Maungawhau 火山、Jura、公的機関の標高データ）は M5-5 で選び、この節に追記する。

#### 8.4.8 評価指標

1. **近似誤差**: RMSE $= \sqrt{\tfrac1M \sum (f_{\mathrm{true}} - f_{\mathrm{interp}})^2}$、最大誤差 $\max \lvert f_{\mathrm{true}} - f_{\mathrm{interp}} \rvert$、収束次数（§8.4.2）。
2. **SciPy との一致**: §8.2 の区分で判定する（一律の許容誤差は使わない）。
3. **スループット**: 構築時間と、評価点数あたりの評価時間を別々に報告する。API は `eval_batch`・`eval_batch_into`・`eval_batch_par`・`eval_batch_into_par`。
4. **メモリ**: 構築時と評価時のピーク確保量を別々に報告する（§8.3）。

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
