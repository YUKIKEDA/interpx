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

* 三重対角系（CubicSpline 等）は **O(n) の Thomas 法を自前実装**する（周期境界は Sherman–Morrison による巡回三重対角）。帯行列系（`BSpline1D`、グリッドの `cubic` / `quintic`）は部分ピボットつきの帯行列 LU を自前実装する（§6.1）。どちらにも `faer` は使わない。
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

  * **crate 内の列ごとの演算**: 公開する演算は上の 3 つ（線形結合に要るもの）だけである。補間器の中には、列ごとの非線形な演算が要るものがある（PCHIP と Akima の傾きは、列ごとの符号の比較や調和平均で決まる。`CubicSpline` の周期性の検査は列ごとの差の大きさを見る）。これらは、封印に使っている crate 内のトレイトに持たせ、公開しない。持つのは「いくつかの値を列ごとに関数で合成する」（2 個のものと、個数が任意のもの。Akima の傾きは 4 本の割線と閾値の 5 個から決まる）と「列ごとの近さの判定」である。利用者が新しい値の型を足すことはないので、公開の `Value` を広げる必要がない。
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

pub enum BoundaryCondition<V> {                            // CubicSpline / BSpline1D。V は値の型
    NotAKnot,                                 // 既定（SciPy に合わせる）
    Periodic,
    Ends(EndCondition<V>, EndCondition<V>),   // 左端、右端
}
pub enum EndCondition<V> { NotAKnot, FirstDerivative(V), SecondDerivative(V) }
```
* **`Extrapolate`**: `Forbid` は領域外で `None`。`Clamp` は端の値を返す（`numpy.interp` の挙動）。`Extend` は端の多項式をそのまま延長する（SciPy の `extrapolate=True`）。`Periodic` は領域に折り返す。
* **既定値は型ごとに SciPy 1.18 に合わせる。**

  | 型                                                                      | 既定        |
  | :---------------------------------------------------------------------- | :---------- |
  | `Linear1D`                                                              | `Clamp`     |
  | `CubicSpline`、`PchipInterpolator`、`CubicHermiteSpline`、`BSpline1D`、`BSplineCurve` | `Extend`    |
  | `CubicSpline`・`BSpline1D` で境界条件が `Periodic` のとき               | `Periodic`  |
  | `Akima1D`、`Nearest1D`、`RegularGridInterpolator`                       | `Forbid`    |

* **`BoundaryCondition`**: 端ごとに条件を持ち、左右で別の条件を指定できる（SciPy の `bc_type=((1, a), (2, b))`）。`Natural` は両端 `SecondDerivative(0)`、`Clamped` は両端 `FirstDerivative(0)` で、`BoundaryCondition::natural()` / `clamped()` で作る。周期は両端まとめての指定で、片側だけ周期という状態は型で表せない。型パラメータは**値の型**である（M1-8 で確定）。端の導関数は値と同じ形を持つので、値が `[T; M]` のスプラインでは `FirstDerivative([..; M])` と列ごとに与える（SciPy も列ごとの配列を受け取る）。
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
* **点ごとの入力が値のほかにもある型**（`CubicHermiteSpline` の傾き）は、上の各段で、それを値の次に検査する。段 1 は値の長さ、次に傾きの長さ。段 3 は座標、値、傾きの順である。段の順は変えない（傾きに NaN があり座標が未ソートなら、`NonFinite` を返す）。

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
* `#[non_exhaustive]` なので、variant は最初に必要になる行で追加する（使う前に定義しない）。上は M1（1D・curve）が使うものである。`NotPeriodic` は M1-2 で定義済みで、最初に返すのは M1-8（`CubicSpline`）である。
* `ParamReason` には、M1-8 が `NonFinite`（引数の値が NaN または無限大。境界条件の導関数）を追加する。
* M1-12 が `TooSmall { min: usize }`（下限より小さい。`BSpline1D` の次数 0）と `Unsupported`（値そのものは正しいが、ほかの引数との組み合わせでは使えない。3 次以外の `BSpline1D` に端の導関数を指定した場合）を追加する。
* 後の行で追加する variant:

  | variant                                                    | 追加する行         | 用途                                                   |
  | :--------------------------------------------------------- | :----------------- | :----------------------------------------------------- |
  | `InvalidParameter { name: &'static str, reason: ParamReason }` | M1-4（`BSpline`） | 引数が範囲外。M1-4 は `ParamReason::TooLarge { max: usize }`（次数の上限）を入れる。必須欠落・最小次数未満などの `ParamReason` は、それを返す行で足す |
  | `GridAxis { axis: usize, fault: AxisFault }` | M2-1（`RegularGridInterpolator`） | グリッドの軸 `axis` の誤り。`AxisFault` は `InsufficientPoints { required, actual }`、`NonFinite { index }`、`NotMonotonic { index }`、`DuplicatePoints { first, second }`（位置はその軸の中のもの）。1D の variant は、どの軸かを持てない |
  | `SolverFailed(SolverError)`                                | M3（RBF）          | 特異な連立方程式（`SolverError::Singular`）            |
  | `DegenerateGeometry(DegeneracyKind)`                       | M4（`LinearND`）   | `LowerDimensional { affine_dim }`、`NumericalLimit`    |

* `NonFinite { index }` の `index` は、座標に非有限値があればその座標の位置で、座標がすべて有限のときに限り値の位置である（検証は座標が先。§3.5）。variant からはどちらの配列かを区別できない。
* 文字列のペイロードは持たない（パターンマッチできる型付きの値だけ）。`ParamReason`・`AxisFault`・`SolverError`・`DegeneracyKind` も `#[non_exhaustive]`。
* 次元の不一致を表す variant は持たない（`const D` のみなので、型で表現できない）。
* `NotSorted { index }` の `index` は、直前の座標より小さい最初の要素。`LengthMismatch` の `expected` は座標（またはクエリ）の長さ、`actual` は値（または出力バッファ）の長さ。
* **`convx::ConvexHullError` は公開しない**（`convx` の版と SemVer を結合させないため）。M4 で、内部で次のように変換する（`From` は `pub(crate)`）:
  * `InsufficientPoints` → `InsufficientPoints`
  * `NonFiniteCoordinate{index}` → `NonFinite{index}`
  * `LengthMismatch` → `LengthMismatch`
  * `DegenerateDimension{actual_dim,..}` → `DegenerateGeometry(LowerDimensional{affine_dim})`
  * `ExactEvaluationExhausted` → `DegenerateGeometry(NumericalLimit)`
  * `TooManyPoints` → `InvalidParameter`
* **確保の大きさ**: 入力だけで決まる上限はエラーで返し、機械のメモリで決まる失敗は捕まえない。前者は、呼び出し側の引数で確保の大きさが決まり、それが `Vec` の持てる上限（`isize::MAX` バイト）を超える場合で、確保の前に検査して `InvalidParameter { reason: TooLarge { max } }` を返す（そのまま確保すると panic するため）。後者（上限内でのメモリ不足）は Rust の標準どおりプロセスが終了する。メモリ不足をエラーとして返す API は持たない。依存クレート（`faer`・`kiddo`・`convx`）の内部の確保は捕まえられず、一部の確保だけをエラーにすると、エラーの値（`max`）が「この値までなら通る」を意味しなくなるためである。
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
│   ├── cubic_hermite # CubicHermiteSpline
│   ├── hermite       # 値と傾きから 3 次の PPoly を作る共通の手順（crate 内。§6.1）
│   ├── bspline1d     # make_interp_spline (次数 1〜15)
│   └── nearest       # Nearest1D (nearest / nearest-up / previous / next)
├── curve/            # N次元パラメトリック曲線
│   └── bspline_curve # BSplineCurve、Parametrization
├── grid/             # N次元直交グリッド補間
│   └── regular_grid  # RegularGridInterpolator
└── scattered/        # N次元不規則散乱データ補間
    ├── nearest_nd    # NearestNDInterpolator (kiddo)
    ├── linear_nd     # LinearNDInterpolator (convx + Jump-and-Walk)
    ├── clough_tocher # CloughTocher2DInterpolator
    └── rbf           # RBFInterpolator (faer + kiddo)
```
（モジュール名は Rust の識別子制約を避け `oned` とする。）

* これはソースの置き場である。`oned` などのモジュールは**非公開**で、公開する項目はすべて `lib.rs` から再輸出する。公開パスは `interpx::Linear1D` のように crate 直下になる（M1-7 で確定）。
* 補間器が共用する基盤（`PPoly`・`BSpline`・区間探索・入力検査・三重対角ソルバー）は、この木の外、`src/` の直下に置く。

### 5.2 網羅対応表

| 分類              | SciPy 1.18                                      | `interpx` 実装構造体                          | 内部コアエンジン                       |
| :---------------- | :---------------------------------------------- | :-------------------------------------------- | :------------------------------------- |
| **1D**            | `numpy.interp`                                  | `Linear1D`                                    | `PPoly`                                |
| **1D**            | `CubicSpline`                                   | `CubicSpline`                                 | Thomas 法 + `PPoly`                    |
| **1D**            | `PchipInterpolator`                             | `PchipInterpolator`                           | `PPoly`（単調性制限導関数）            |
| **1D**            | `Akima1DInterpolator`（`method='makima'` 含む） | `Akima1D`                                     | `PPoly`                                |
| **1D**            | `CubicHermiteSpline`                            | `CubicHermiteSpline`                          | `PPoly`                                |
| **1D**            | `make_interp_spline`（次数 1〜15）              | `BSpline1D`                                   | `BSpline`（de Boor）                   |
| **1D**            | `interp1d` の nearest / previous / next         | `Nearest1D`                                   | 二分探索                               |
| **1D**            | `PPoly` / `BSpline` の導関数・積分              | `derivative` / `antiderivative` / `integrate` | `PPoly` ⇄ `BSpline` 相互変換           |
| **N-D Curve**     | `make_interp_spline`（N-D y）                   | `BSplineCurve<T, D>`                          | `BSpline`（ベクトル制御点）            |
| **N-D Grid**      | `RegularGridInterpolator` `nearest`             | `RegularGridInterpolator`                     | インデックス四捨五入                   |
| **N-D Grid**      | `linear` / `slinear`（同じ関数）                | 同上（`GridMethod::Linear`）                  | テンソル積線形                         |
| **N-D Grid**      | `cubic` / `quintic`                             | 同上                                          | 構築時テンソル積 B-spline 係数（§6.4） |
| **N-D Grid**      | `pchip`                                         | 同上                                          | 軸ごとの PCHIP を最後の軸から（§6.4） |
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
  * **領域外**: `Extrapolate` の 4 モード。既定は `Extend`（SciPy の `PPoly` と同じ）で、`with_extrapolate` で変える。`Periodic` は $x_0 + ((x - x_0) \bmod (x_n - x_0))$ に折り返す。右端 $x_n$ も $x_0$ に折り返す（1 周期は左端だけを含む。SciPy と同じ）。値が周期的な関数では差が出ないが、両端で値や導関数が違うときは左端の区間の値になる。NaN のクエリは `None`。無限大のクエリは `Forbid` と `Periodic` で `None`、`Clamp` は端の値、`Extend` は端の多項式をそのまま評価した値（一般に無限大か NaN）。
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
  * **次数の上限**: 評価は確保しない（§3.2）。作業用の配列をスタックに固定長で置くため、次数は 15 までとする。SciPy の `make_interp_spline` が実用上使うのは 5 次までである。上限はクレート直下の定数 `BSPLINE_MAX_DEGREE` で公開する（`BSpline` の関連定数にすると、型パラメータを指定しないと書けない）。
  * **領域外**: `Extrapolate` の 4 モードで、既定は `Extend`（SciPy の `BSpline` と同じ）。`Extend` は最初または最後の幅が正の区間の多項式を延長する。`Periodic`・NaN・無限大のクエリの扱いは `PPoly` と同じ。
  * **区間探索の共有**: 二分探索、`Hint` つきの探索、バッチで区切りごとにヒントを使うかの判定は、crate 内のモジュール `search` に置き、`PPoly` と `BSpline` が共用する。`Hint` 型も共通で、`BSpline` も `eval_hinted` を持つ。
* **`PPoly` と `BSpline` の変換**: どちらも厳密な変換（丸め誤差を除く）で、外挿モードを引き継ぐ。
  * **`PPoly::from_bspline(&BSpline) -> Result<PPoly, _>`**（SciPy の `PPoly.from_spline`）: 区切り点は定義域 $[t_k, t_n]$ 内の相異なるノット。区間 $[u_j, u_{j+1}]$ の係数は、スプラインの $m$ 階導関数の $u_j$ での値を $m!$ で割ったもの（$u_j$ まわりのテイラー係数）。SciPy は全ノットを区切り点に残すので重複ノットが幅ゼロの区間になるが、`PPoly` の区切り点は狭義単調増加なので、幅ゼロの区間は作らない。
  * **`BSpline::from_ppoly(&PPoly) -> Result<BSpline, _>`**: 各区切り点に多重度 $k+1$ のノットを置き、区間ごとにべき基底の係数を Bernstein 基底の係数（Bézier の制御点）に直す。区間どうしの連続性を仮定しないので、不連続な区分多項式でも厳密に表せる。制御点は `区間数 × (k + 1)` 個になる。連続性を調べてノットの多重度を減らす処理は v0.1.0 の対象外（ロードマップ A-9）。次数が 15 を超える `PPoly` は `InvalidParameter` を返す。**領域外での精度**: 区切り点の内側では 2 つの表現は丸めの範囲で一致する。端の区間が短いとき、その Bézier 表現を領域の外へ遠く延ばすと、`PPoly` が保つ桁を失う（端の区間の幅の $r$ 倍先で、丸め誤差はおよそ $r^k$ 倍になる。幅 $10^{-3}$ の 3 次の区間を 700 幅先で評価して相対 $4 \times 10^{-8}$）。遠い外挿は `PPoly` で評価する。
  * 往復（`PPoly` → `BSpline` → `PPoly`、およびその逆）で、評価値が許容誤差内で一致することをテストする。
* **導関数と積分（M1-5）**: `PPoly` と `BSpline` に、同じ形の 4 つの関数を置く。SciPy 1.18 は 3 か所で 2 つの型の挙動が違うが、ここでは統一する（下の「SciPy との差」）。
  * **関数**:
    * `derivative(nu: usize) -> Result<Self, InterpolateError>`: $\nu$ 階導関数。
    * `antiderivative(nu: usize) -> Result<Self, InterpolateError>`: $\nu$ 階の原始関数。
    * `integrate(a: T, b: T) -> Option<V>`: $a$ から $b$ までの定積分。
    * `eval_derivative(x: T, nu: usize) -> Option<V>`: $x$ での $\nu$ 階導関数の値。
  * `derivative` と `antiderivative` は新しいオブジェクトを作る（確保する）。`nu = 0` は複製を返す。係数があふれたら `NonFinite` を返す。`integrate` と `eval_derivative` は評価であり、確保しない（§3.2。`tests/alloc.rs` の対象）。
  * **`derivative`**:
    * `PPoly`: 区切り点はそのまま、次数は $k - \nu$。$(x - x_i)^m$ の係数に $m(m-1)\cdots(m-\nu+1)$ を掛け、低次の $\nu$ 個を落とす。
    * `BSpline`: 1 階ごとに、ノットの両端を 1 個ずつ落とし、制御点を $c'_i = k\,\dfrac{c_{i+1} - c_i}{t_{i+k+1} - t_{i+1}}$（$n - 1$ 個）に置き換える。定義域 $[t_k, t_n]$ は変わらない。分母が 0 になるのは $t_{i+1} = \dots = t_{i+k+1}$ のときで、その基底関数は台が空なので、制御点は 0 とする。これで、`from_ppoly` が作る多重度 $k+1$ のノットを持つスプラインも微分できる。
    * $\nu > k$ のときは、どちらの型も恒等的に 0 を返す。`PPoly` は区切り点をそのまま持つ 0 次、`BSpline` はノットが定義域の両端だけの 0 次である。
  * **`antiderivative`**:
    * `PPoly`: 次数は $k + \nu$。1 階ごとに、$(x - x_i)^m$ の係数を $m + 1$ で割って次数を上げ、定数項は、原始関数が区切り点で連続になるように左から足し上げる。`PPoly` に次数の上限はないが、結果の係数が `Vec` の持てる上限（`isize::MAX` バイト）を超える `nu` は `InvalidParameter { name: "nu", reason: TooLarge { max } }` を返す（`max` はその上限に収まる最大の階数。確保の扱いは §4）。
    * `BSpline`: 1 階ごとに、ノットの両端に端のノットを 1 個ずつ足し、制御点を $C_0 = 0$、$C_{j+1} = C_j + c_j\,\dfrac{t_{j+k+1} - t_j}{k+1}$（$n + 1$ 個）とする。結果の次数が 15 を超えるときは `InvalidParameter { name: "nu", reason: TooLarge { max } }` を返す（`max` は $15 - k$）。
    * **積分定数**: どちらの型も、定義域の左端で 0 にする。`BSpline` は上の制御点から、定義域の左端での値を引く（基底関数の和が 1 なので、全制御点から同じ値を引けばよい）。端のノットが $k+1$ 重なら引く値は厳密に 0 で、SciPy と一致する。
  * **引き継ぐ外挿モード**:

    | 元のモード | `derivative` | `antiderivative` |
    | :--------- | :----------- | :--------------- |
    | `Extend`   | `Extend`     | `Extend`         |
    | `Forbid`   | `Forbid`     | `Forbid`         |
    | `Periodic` | `Periodic`   | `Forbid`         |
    | `Clamp`    | `Forbid`     | `Forbid`         |

    周期関数の原始関数は一般に周期的でない（SciPy も「外挿なし」に変える）。`Clamp` は領域外を定数で延ばすので、導関数は外側で 0、原始関数は外側で 1 次式になり、どちらも 4 つのモードでは表せない。表せない領域は値を返さない `Forbid` にする。`Clamp` の外側の値が要るときは、元のオブジェクトの `eval_derivative` と `integrate` が返す。
  * **`eval_derivative`**: `eval` が定める関数の $\nu$ 階導関数を返す。
    * 区間の決め方は `eval` と同じ。導関数が跳ぶ区切り点では右側の区間の値、定義域の右端だけは左側の区間の値になる（SciPy と同じ）。
    * 領域外は、`Extend` が端の多項式の導関数、`Periodic` が折り返した点での値、`Forbid` が `None`。`Clamp` は $\nu \ge 1$ なら 0（そこでは定数）。NaN は `None`。
    * $\nu = 0$ は `eval` とビット同一。$\nu > k$ は、値を持つ点で 0。
    * 評価の式は、`PPoly` が係数に階乗の比を掛けながらのホーナー法、`BSpline` が導関数用の漸化式（上の「漸化式の形」）。
    * **バッチ版とヒント版は置かない**。多数の点では `derivative(nu)` で作ったオブジェクトに `eval_batch_into` を使う。`BSpline` はこの経路のほうが 1 点あたりも速い（値の漸化式のほうが速いため）。`derivative(nu)` は `Clamp` を引き継がないので、`Clamp` の外側の点は `eval_derivative` で求める。
  * **`integrate`**: `eval` が定める関数の定積分を返す。
    * $a > b$ は符号を反転した値、$a = b$ は 0。$a$ か $b$ が有限でなければ `None`。
    * `Forbid`: $a$ か $b$ が定義域の外なら `None`。
    * `Extend`: 端の多項式を延ばして積分する。
    * `Clamp`: 外側の部分は、端の値 × 長さ。
    * `Periodic`: 周期を $P$ として、$b - a$ を $P$ の整数倍と残りに分け、整数倍の分は 1 周期の積分の倍数、残りは折り返した区間（定義域の右端をまたぐときは 2 つに分ける）で求める（SciPy と同じ手順）。
    * `PPoly`: またぐ区間ごとに、多項式の原始関数の両端での差を足す。計算量はまたぐ区間の数に比例する。
    * `BSpline`: 原始関数 $F$ の制御点を配列に作らず、$F(b) - F(a)$ を直接求める。$a$ の区間 $\ell_a$ と $b$ の区間 $\ell_b$ で、次数 $k+1$ の基底関数の値（各 $k+2$ 個）をスタックに求め、$j = \ell_a - k$ から $\ell_b + 1$ まで $C_j$ を足し上げながら、重み（$b$ での基底の値 − $a$ での基底の値）を掛けて足す。基底関数の和は 1 なので、$C_{\ell_a - k} = 0$ から始めても差は変わらない。計算量はまたぐノット区間の数に比例する。次数 $k+1$ の基底は元のノット列の上でそのまま計算できる（足した両端のノットは、定義域の区間の漸化式に現れない）。作業配列は 17 個で、次数 15 のスプラインも積分できる。
  * **SciPy との差**（`PPoly` と `BSpline` で挙動をそろえるため）:

    | 場面 | SciPy 1.18 | `interpx` |
    | :--- | :--------- | :-------- |
    | `BSpline` の $\nu > k$ の `derivative` | `ValueError` | 恒等的に 0（`PPoly` と同じ） |
    | 内部に多重度 $k+1$ のノットがある `BSpline` の `derivative` | `ValueError` | 微分できる |
    | 端のノットが $k+1$ 重でない `BSpline` の原始関数 | 定義域の左端で 0 にならない | 定義域の左端で 0（`PPoly` と同じ） |
    | 外挿なしで定義域を出る `integrate` | `PPoly` は NaN、`BSpline` は定義域に切り詰めて積分 | どちらも `None` |
    | NaN の積分端 | `PPoly` は例外、`BSpline` は NaN | どちらも `None` |

    いずれも SciPy 1.18.1 を実行して確かめた。
* **三重対角ソルバー（M1-6）**: `CubicSpline` の構築が解く連立方程式のためのもので、crate 内に置き公開しない（`scipy.interpolate` にも対応する公開関数はない）。
  * **行列の持ち方**: 長さ $n$ の 3 本の配列 `lower`・`diag`・`upper` で、`lower[i]` が $A_{i,i-1}$、`diag[i]` が $A_{i,i}$、`upper[i]` が $A_{i,i+1}$。`lower[0]` と `upper[n-1]` は三重対角の外で、非周期版は読まない。
  * **右辺**: `V: Value<T>` の配列。値が `[T; M]` のとき、$M$ 列を 1 回の消去で解く。スカラーは `T: Scalar`。
  * **確保しない**: 渡された配列を書き換える。`diag` は消去後のピボットの逆数に、右辺は解に置き換わる。`lower` と `upper` は読むだけである。構築時の確保は呼び出し側の 1 か所にまとまる。
  * **解法**: ピボット選択のない Thomas 法（前進消去と後退代入、$O(n)$）。下りの 1 回で消去と前進代入を、上りの 1 回で後退代入を行う。巡回版は右辺が 2 本（データの右辺と、Sherman–Morrison の $u$）で、同じ 2 回の走査に両方を乗せる。
  * **依存の鎖を短くする**: この解法は漸化式で、各行が直前の行の結果を待つ。速度を決めるのは演算の個数より、ある行の結果から次の行の結果までに並ぶ演算の列（依存の鎖）の長さである。ピボットは $p_i = d_i - (l_i u_{i-1})/p_{i-1}$ と書き、積 $l_i u_{i-1}$ を鎖の外で求めて、鎖を除算 1 回と減算 1 回にする。代入は、ピボットの逆数を先に掛けてから隣の行を引く形にして、鎖を乗算 1 回と加算 1 回にする。消去と代入を別々の走査にすると鎖の時間が足し合わされるので、1 回の走査にまとめる。「除算を 1 行 1 回に減らす」形（ピボットの逆数だけを漸化式で持つ）は、除算は減るが鎖が長くなり、速くならなかった（`benches/cubic_spline.rs` の実測。1000 点の構築で 21 µs のまま。鎖を短くした形は 16 µs）。
  * **ピボット選択をしない理由**: SciPy の `CubicSpline` は部分ピボットつきの帯行列ソルバー（`solve_banded`）で解くが、ここで解く行列は選択なしで安定に解ける。自然・固定（1 階導関数）・2 階導関数の境界条件の行列は狭義優対角である。not-a-knot は最初と最後の行が優対角でないが、最初の行 $[\,h_1,\; h_0 + h_1\,]$ で 2 行目 $[\,h_1,\; 2(h_0 + h_1),\; h_0\,]$ を消去したあとのピボットは $h_0 + h_1 > 0$ である（$h_i = x_{i+1} - x_i$）。最後の行 $[\,h_{n-3} + h_{n-2},\; h_{n-3}\,]$ のピボットは $h_{n-3}\,\bigl(1 - (h_{n-3} + h_{n-2})/p\bigr)$ で、直前の行のピボット $p$ は $2h_{n-3} + h_{n-2}$ より大きいので、これも正である。SciPy との一致は、M1-8 で `CubicSpline` の golden にかけて確かめる。
  * **破綻**: ピボットが 0 または非有限になったら、その行を持つ crate 内のエラーを返す（panic しない）。公開のどの variant に写すかは、最初にそれを返す M1-8 で決める。解が有限かどうかは呼び出し側が検査する（`PPoly` の構築が係数を検査する）。
  * **巡回三重対角（周期境界）**: 三重対角に加えて右上 $A_{0,n-1}$ と左下 $A_{n-1,0}$ に 1 要素ずつを持つ系を、Sherman–Morrison の公式で解く。右上は `lower[0]`、左下は `upper[n-1]` に置く（非周期版が読まない 2 要素を使う）。$\gamma = -A_{0,0}$ として、$A = A' + u v^{\mathsf T}$、$u = (\gamma, 0, \dots, 0, A_{n-1,0})^{\mathsf T}$、$v = (1, 0, \dots, 0, A_{0,n-1}/\gamma)^{\mathsf T}$ と分ける。$A'$ は三重対角で、$A$ の対角の両端を $A_{0,0} - \gamma$ と $A_{n-1,n-1} - A_{n-1,0}A_{0,n-1}/\gamma$ に替えたものである。$A' y = b$ と $A' z = u$ を解き、$x = y - z\,\dfrac{v^{\mathsf T} y}{1 + v^{\mathsf T} z}$ とする。消去は 1 回、代入は 2 回で、$z$ 用の長さ $n$ の作業配列は呼び出し側が渡す。
  * **大きさ**: 非周期版は $n \ge 1$。巡回版は $n \ge 3$（$n \le 2$ では巡回の形が定まらない）。点数が少ない周期スプラインは、SciPy と同じく `CubicSpline` 側で別に扱う（M1-8）。
* **帯行列ソルバー（M1-12）**: `BSpline1D` の選点行列と、グリッドの `cubic` / `quintic`（§6.4）の軸ごとの系は、三重対角より帯が広い（次数 $k$ で、not-a-knot のノット列では下帯・上帯とも最大 $k$）。これは**部分ピボットつきの帯行列 LU** で解く（SciPy は LAPACK の `gbsv`）。M1-12 が crate 内に実装し、M2-2 がそれを使う。三重対角を帯行列 LU で兼ねないのは、`CubicSpline` の構築がもっともよく使われる経路で、専用の Thomas 法のほうが速いためである。
  * **場所と形**: crate 内のモジュール `banded`。行列は `T`、右辺は `V: Value<T>` なので、値が `[T; M]` でも列ごとに独立に解ける（ピボットの選択は行列だけで決まる）。分解と求解を分ける: 分解は 1 回で、右辺は型の違うものを何本でも解ける（周期スプラインが、`T` の右辺と `V` の右辺を同じ分解で解く）。
  * **格納**: 下帯 $k_l$ 本、上帯 $k_u$ 本の $n \times n$ 行列を、LAPACK の `gbsv` と同じ帯形式で持つ（1 列あたり $2k_l + k_u + 1$ 個。上の $k_l$ 本は、行の入れ替えで生じる埋まりのための空きである）。$k_l = k_u = n - 1$ とすれば密行列も同じ関数で解ける（周期スプラインの小さな系に使う）。
  * **消去**: 列ごとに、対角以下で絶対値が最大の要素をピボットに選び、行を入れ替える（部分ピボット）。ピボットが 0 か非有限なら `Breakdown { row }` を返す（三重対角ソルバーと同じ形。panic しない）。解が有限かどうかは呼び出し側が確かめる。
  * **確保**: 行列と入れ替えの記録は呼び出し側が渡した領域に置く。ソルバーの中では確保しない。
  * **巡回帯行列**: 帯の外に、右上と左下の隅のブロックを持つ行列（周期スプライン）は、隅を除いた帯行列 $A$ と低ランクの補正 $U V^{\mathsf T}$ に分け、Woodbury の公式で解く（SciPy の `_woodbury_algorithm` と同じ）: $A Z = U$ と $A y = b$ を解き、$x = y - Z\,(I + V^{\mathsf T} Z)^{-1} V^{\mathsf T} y$。補正の階数は $k_l + k_u$ で、$I + V^{\mathsf T} Z$ はその大きさの密行列である（同じソルバーで解く）。隅が帯と重ならないこと（$n > k_l + k_u$）が前提で、それより小さい系は密行列として解く。
* **区分多項式である補間器と `PPoly`**: `Linear1D`・`CubicSpline`・`PchipInterpolator`・`Akima1D`・`CubicHermiteSpline` は、内部に `PPoly` を持ち、`as_ppoly(&self) -> &PPoly<T, V>` で貸し出す（確保しない）。導関数と積分（`derivative`・`antiderivative`・`integrate`・`eval_derivative`）は、この `PPoly` のものを使う。補間器の型に同じ 4 つの関数を重ねて持たせない。`PPoly` の外挿モードは補間器のものと同じにする。 `BSpline1D`（M1-12）は同じ規則で `BSpline` を持ち、`as_bspline(&self) -> &BSpline<T, V>` で貸し出す。
* **`Linear1D<T: Scalar, V: Value<T> = T>`（M1-7、`numpy.interp`）**: 折れ線。
  * **構築**: `Linear1D::new(&[T], &[V])`。検査は §3.5 の順で、点は 2 個以上（`InsufficientPoints { required: 2 }`）。
  * **保持するもの**: 1 次の `PPoly`（区間 $i$ の係数は、傾き $s_i = (y_{i+1} - y_i)/(x_{i+1} - x_i)$ と左端の値 $y_i$）と、最後のデータ値 $y_{n-1}$。
  * **評価**: 区間 $i$ で $y_i + s_i\,(x - x_i)$（`numpy.interp` と同じ式）。左端 $x = x_i$ では $y_i$ が厳密に返る。区間の右端 $x = x_{i+1}$ に当たるのは定義域の右端だけで、そこでは式を使わず $y_{n-1}$ をそのまま返す（式のままだと 1 ulp ずれることがある）。したがって、すべての節点でデータ値が厳密に返る。`as_ppoly()` の `eval` は式のままなので、定義域の右端で 1 ulp 違うことがある。
  * **領域外**: 4 モード。既定は `Clamp`（`numpy.interp` と同じく端のデータ値を返す。右外は $y_{n-1}$ が厳密に返る）。`Extend` は端の線分を延ばす。`Periodic` は他の型と同じ、定義域への折り返しである。
  * **`numpy.interp` の引数との対応**: 外側の埋め値 `left` / `right` は持たない。`Forbid` にして `eval_or(x, fill)` を使う（左右で別の値にするなら、呼び出し側が `None` を受けて定義域と比べる）。`period` も持たない。`numpy.interp` の `period` は節点も周期で折り返して並べ直すもので、`Periodic` とは意味が違う。
  * **`eval_hinted` とバッチ**: `PPoly` と同じ（`search` を使う）。
* **`Nearest1D<T: Scalar, V: Value<T> = T>`（M1-7、`interp1d` の `nearest` / `nearest-up` / `previous` / `next`）**: 最も近い（または直前・直後の）データ点の値を返す。
  * **構築**: `Nearest1D::new(&[T], &[V], NearestKind)`。`NearestKind::{Nearest, NearestUp, Previous, Next}`（`#[non_exhaustive]`）。検査は §3.5 の順で、点は 2 個以上。
  * **保持するもの**: 座標と値と種類だけ。`PPoly` は持たず、`as_ppoly()` もない（区分定数だが、点がどの区間に属するかが `PPoly` と違う）。
  * **どの点の値を返すか**: 定義域内の $x$ が区間 $[x_i, x_{i+1}]$ にあるとき、次の点の値を返す。

    | 種類 | 返す点 |
    | :--- | :----- |
    | `Nearest` | 中点 $m_i$ に対し、$x \le m_i$ なら $x_i$、そうでなければ $x_{i+1}$（中点ちょうどは左） |
    | `NearestUp` | $x < m_i$ なら $x_i$、そうでなければ $x_{i+1}$（中点ちょうどは右） |
    | `Previous` | $x < x_{i+1}$ なら $x_i$、$x = x_{i+1}$ なら $x_{i+1}$（$x$ 以下で最大の節点） |
    | `Next` | $x = x_i$ なら $x_i$、そうでなければ $x_{i+1}$（$x$ 以上で最小の節点） |

    中点は SciPy と同じ式 $m_i = x_i/2 + x_{i+1}/2$ で求める（和を先に取らないので、あふれない）。丸めまで SciPy と同じ点を選ぶ。
  * **領域外**: 4 モード。既定は `Forbid`（SciPy の `interp1d` は既定で例外を出す）。`Clamp` と `Extend` はどちらも端のデータ値を返す（区分定数を延ばしても同じ値）。`Periodic` は折り返す。NaN のクエリは `None`。
  * **`eval_hinted` とバッチ**: `search` で区間を求め、上の表で点を選ぶ。`PPoly` と同じ経路なので、ヒントとバッチの結果は `eval` とビット同一である。
* **NumPy・SciPy との差（M1-7）**:

  | 場面 | NumPy 2.5 / SciPy 1.18 | `interpx` |
  | :--- | :--------------------- | :-------- |
  | データ点が 1 個 | `numpy.interp` は定数、`interp1d` も構築できる | `InsufficientPoints`（区間探索は 2 点以上を前提にする） |
  | 重複した $x$、並んでいない $x$ | `numpy.interp` は検査しない。`interp1d` は並べ替え、重複を受け付ける | `DuplicatePoints` / `NotSorted`（§3.5） |
  | 値の無限大 | `numpy.interp` は受け付ける | `NonFinite`（§3.5） |
  | `Nearest1D` の NaN のクエリ | 最後のデータ値 | `None` |
  | `fill_value="extrapolate"` の `previous` の左外、`next` の右外 | NaN | `Clamp` / `Extend` は端のデータ値 |

  いずれも NumPy 2.5.4 と SciPy 1.18.1 を実行して確かめた。
* **`CubicSpline<T: Scalar, V: Value<T> = T>`（M1-8、SciPy の `CubicSpline`）**: $C^2$ の 3 次スプライン。
  * **構築**: `CubicSpline::new(&[T], &[V], BoundaryCondition<V>)`。境界条件は必須の引数で、既定の not-a-knot は `BoundaryCondition::default()` で書く。点は 2 個以上。
  * **構築時の確保**: 3 個（座標の複製、係数、作業領域）。傾きは係数の配列に解き、係数を最後の区間から前へ向かって上書きする（区間 $i$ は傾き $s_i$・$s_{i+1}$ を読んで位置 $4i$ 以降に書くので、まだ要る傾きを壊さない）。作業領域は、行列の 3 本の対角と、区間の幅の逆数（割線と係数の両方で使うので 1 回だけ求める）、巡回版の作業配列を 1 つの配列に持つ。入力検査は 1 回で、検査済みの部品から `PPoly` を作る crate 内の構築経路（`PPoly::from_checked`）を使う。`Linear1D` も同じ経路を使う。
  * **保持するもの**: 3 次の `PPoly`。`as_ppoly()` で貸し出す。評価・`eval_hinted`・バッチは `PPoly` のものをそのまま使う（`Linear1D` と違い、右端の特別扱いはしない。SciPy も右端は多項式の値を返す）。
  * **求め方**: 各データ点での傾き $s_i$ を連立方程式で求め、区間 $i$ の係数を傾きから作る。$h_i = x_{i+1} - x_i$、$\delta_i = (y_{i+1} - y_i)/h_i$ とおく。
    * 内側の行（$1 \le i \le n-2$）: $h_i\,s_{i-1} + 2(h_{i-1} + h_i)\,s_i + h_{i-1}\,s_{i+1} = 3\,(h_i\,\delta_{i-1} + h_{i-1}\,\delta_i)$。
    * 端の行（左端。右端は対称）:

      | 端の条件 | 行 | 右辺 |
      | :------- | :- | :--- |
      | `NotAKnot` | $h_1\,s_0 + (h_0 + h_1)\,s_1$ | $\bigl((h_0 + 2d)\,h_1\,\delta_0 + h_0^2\,\delta_1\bigr)/d$、$d = h_0 + h_1$ |
      | `FirstDerivative(v)` | $s_0$ | $v$ |
      | `SecondDerivative(v)` | $2h_0\,s_0 + h_0\,s_1$ | $-\tfrac12 v\,h_0^2 + 3\,(y_1 - y_0)$ |

    * 係数: $t = (s_i + s_{i+1} - 2\delta_i)/h_i$ として、区間 $i$ は $\dfrac{t}{h_i}(x - x_i)^3 + \Bigl(\dfrac{\delta_i - s_i}{h_i} - t\Bigr)(x - x_i)^2 + s_i\,(x - x_i) + y_i$。
    * 系は三重対角ソルバー（Thomas 法、ピボット選択なし）で解く。行列と右辺の組み立ては SciPy（`_cubic.py`）と同じである。
  * **周期境界**: $s_{n-1} = s_0$ なので、未知数は $s_0, \dots, s_{n-2}$ の $n-1$ 個。行 $i$ は内側の行と同じ式で、添字を $n-1$ を法として回す。右上と左下に 1 要素ずつが出る巡回三重対角になり、巡回ソルバー（Sherman–Morrison）で解く。
  * **点が少ない場合**（SciPy と同じ扱い）:

    | 点数 | 境界条件 | 結果 |
    | :--- | :------- | :--- |
    | 2 | `NotAKnot`（どちらの端でも） | その端の条件を `FirstDerivative(δ_0)` に置き換える。両端なら直線 |
    | 2 | `Periodic` | 両端の条件を `FirstDerivative(δ_0)` に置き換える。$y_0 = y_1$ なので、定数になる |
    | 3 | 両端 `NotAKnot` | 3 点を通る放物線。行は $[1, 1, 0]$、$[h_1, 2(h_0 + h_1), h_0]$、$[0, 1, 1]$、右辺は $2\delta_0$、$3(h_0\delta_1 + h_1\delta_0)$、$2\delta_1$（三重対角なので同じソルバーで解く） |
    | 3 | `Periodic` | 未知数が 2 個で、巡回ソルバーの範囲外。$s_0 = s_1 = s_2 = \dfrac{\delta_0/h_0 + \delta_1/h_1}{1/h_0 + 1/h_1}$ |
    | 4 以上 | `Periodic` | 巡回ソルバー（3 行以上） |

    上以外（3 点で片端だけ `NotAKnot` など）は一般の式のままである。
  * **周期性の検査**: `Periodic` では $y_0$ と $y_{n-1}$ が一致しなければ `NotPeriodic` を返す。一致の判定は SciPy と同じ許容誤差で、列ごとに $|y_0 - y_{n-1}| \le \tau\,(1 + |y_{n-1}|)$、$\tau = 4.5\,\varepsilon$（`f64` で約 $10^{-15}$。SciPy は `rtol = atol = 1e-15`）。厳密な一致にしないのは、`sin(0)` と `sin(2π)` のような丸めの差で拒否しないためである。検査を通ったあとは、$y_{n-1}$ も与えられた値のまま使う（SciPy と同じ。最後の区間の $\delta_{n-2}$ に入る）。
  * **領域外**: 既定は `Extend`、境界条件が `Periodic` なら `Periodic`（§3.3、SciPy と同じ）。`with_extrapolate` で 4 モードのどれにも変えられる。
  * **エラー**: 検査の順は、§3.5 の共通の検査 → 境界条件の導関数の値（NaN・無限大なら `InvalidParameter { name: "bc", reason: ParamReason::NonFinite }`）→ 周期性（`NotPeriodic`）→ 求解。三重対角ソルバーの破綻と、傾き・係数の桁あふれは `NonFinite { index }` で返す（`index` は該当するデータ点）。有効な入力に対して行列は正則なので、破綻は桁あふれでしか起きない。専用の variant は足さない（`Linear1D` の傾きの桁あふれと同じ扱い）。
  * **SciPy との差**: 共通の検査による差（重複・未ソートの $x$、非有限の値）だけである。点数の下限（2 点）、少ない点の扱い、周期性の許容誤差は SciPy と同じ。
* **値と傾きから 3 次の `PPoly` を作る（3 次エルミート形。M1-9 で切り出す）**: `CubicSpline`・`PchipInterpolator`・`Akima1D`・`CubicHermiteSpline` は、各データ点の値 $y_i$ と傾き $s_i$ が決まれば、同じ手順で 3 次の `PPoly` になる。違うのは傾きの決め方だけである。この手順を crate 内の関数 1 つにし、4 つの型が共用する。
  * 区間 $i$ の係数は、$\delta_i = (y_{i+1} - y_i)/h_i$、$t = (s_i + s_{i+1} - 2\delta_i)/h_i$ として、$\dfrac{t}{h_i}(x - x_i)^3 + \Bigl(\dfrac{\delta_i - s_i}{h_i} - t\Bigr)(x - x_i)^2 + s_i\,(x - x_i) + y_i$（SciPy の `CubicHermiteSpline` と同じ式）。
  * 係数の配列（区間あたり 4 個）の先頭に傾きを置いて呼び、最後の区間から前へ向かって係数で上書きする（M1-8 で `CubicSpline` に入れた形）。確保は増えない。入力検査は呼び出し側が済ませてあり、`PPoly::from_checked` で作る。
  * 係数が有限でなければ `NonFinite { index }`（`index` は最初にあふれた区間）。
* **`PchipInterpolator<T: Scalar, V: Value<T> = T>`（M1-9、SciPy の `PchipInterpolator`）**: 形を保つ $C^1$ の 3 次補間。データが単調な区間では補間も単調で、データの値の範囲を行き過ぎない。
  * **構築**: `PchipInterpolator::new(&[T], &[V])`。検査は §3.5 の順で、点は 2 個以上。
  * **保持するもの**: 3 次の `PPoly`。`as_ppoly()` で貸し出す。評価・`eval_hinted`・バッチは `PPoly` のものである。
  * **傾き**（SciPy と同じ。$h_i = x_{i+1} - x_i$、$m_i = (y_{i+1} - y_i)/h_i$）:
    * **内側の点** $x_i$: $m_{i-1}$ と $m_i$ の符号が違うか、どちらかが 0 なら $s_i = 0$。そうでなければ、重みつきの調和平均 $\dfrac{1}{s_i} = \dfrac{1}{w_1 + w_2}\Bigl(\dfrac{w_1}{m_{i-1}} + \dfrac{w_2}{m_i}\Bigr)$、$w_1 = 2h_i + h_{i-1}$、$w_2 = h_i + 2h_{i-1}$（Fritsch–Butland）。
    * **端**（左端。右端は対称）: 3 点の片側の式 $d = \dfrac{(2h_0 + h_1)\,m_0 - h_0\,m_1}{h_0 + h_1}$ を求め、形を保つように制限する。$d$ の符号が $m_0$ と違えば $s_0 = 0$。そうでなく、$m_0$ と $m_1$ の符号が違い、かつ $|d| > 3|m_0|$ なら $s_0 = 3m_0$。それ以外は $s_0 = d$（Moler の `pchiptx`）。
    * **2 点**: $s_0 = s_1 = m_0$（直線）。
    * 符号は SciPy の `sign` と同じく、正・負・0 の 3 値で比べる。
  * **値が `[T; M]` のとき**: 上の規則を列ごとに独立に当てる（列ごとに符号が違いうる）。§3.1 の「crate 内の列ごとの演算」を使う。
  * **領域外**: 既定は `Extend`（§3.3、SciPy と同じ）。`with_extrapolate` で 4 モードのどれにも変えられる。単調性と「行き過ぎない」性質は領域内のもので、端の 3 次式を延ばした先には及ばない。
  * **エラー**: §3.5 の共通の検査。割線・傾き・係数の桁あふれは `NonFinite { index }`。
  * **SciPy との差**: 共通の検査による差（重複・未ソートの $x$、非有限の値）だけである。
* **`Akima1D<T: Scalar, V: Value<T> = T>`（M1-10、SciPy の `Akima1DInterpolator`）**: $C^1$ の 3 次補間。各点の傾きを近くの 5 点だけから決めるので、外れた 1 点の影響が遠くへ及ばず、スプラインのような波打ちが出にくい。
  * **構築**: `Akima1D::new(&[T], &[V], AkimaMethod)`。`AkimaMethod::{Akima, Makima}`（`#[non_exhaustive]`。SciPy の `method="akima"` / `"makima"`）。検査は §3.5 の順で、点は 2 個以上。
  * **保持するもの**: 3 次の `PPoly`。`as_ppoly()` で貸し出す。係数は共通の手順（3 次エルミート形）で作る。
  * **傾き**（SciPy と同じ）: 割線 $m_i = (y_{i+1} - y_i)/h_i$（$i = 0, \dots, n-2$）を、両端の外へ 2 本ずつ直線で延ばす: $m_{-1} = 2m_0 - m_1$、$m_{-2} = 2m_{-1} - m_0$、$m_{n-1} = 2m_{n-2} - m_{n-3}$、$m_n = 2m_{n-1} - m_{n-2}$。点 $x_i$ の傾きは、その前の 2 本 $m_{i-2}, m_{i-1}$ と後の 2 本 $m_i, m_{i+1}$ から決める。
    * 重み: $f_1 = w(m_i, m_{i+1})$、$f_2 = w(m_{i-2}, m_{i-1})$。`Akima` は $w(a, b) = |b - a|$、`Makima` は $w(a, b) = |b - a| + \tfrac12 |a + b|$。
    * $f_1 + f_2$ が切り捨ての閾値を超える点: $s_i = m_{i-1} + \dfrac{f_2}{f_1 + f_2}\,(m_i - m_{i-1})$。
    * 超えない点（前の 2 本が等しく、後の 2 本も等しい点。傾きが式で定まらない）: $s_i = \tfrac12\,(m_{i-2} + m_{i+1})$。
    * 閾値: $f_1 + f_2$ の全点での最大値の $10^{-9}$ 倍（SciPy の `break_mult`）。
    * **2 点**: $s_0 = s_1 = m_0$（直線）。
  * **`Makima` の性質**: 平らな区間（等しい値が 3 点以上続く）の両端で傾きが 0 になり、平らな区間が平らなまま保たれる。`Akima` は 2 つの割線が等しい点の扱いに飛びがあり（閾値はその位置を動かすだけで、なくせない）、平らな区間の隣で行き過ぎることがある。
  * **値が `[T; M]` のとき**: 傾きも閾値も列ごとに決め、各列を独立に補間する。1 列ずつ作った結果とビット同一になる。
  * **領域外**: 既定は `Forbid`（§3.3。SciPy の Akima は既定で外挿しない）。`with_extrapolate` で 4 モードのどれにも変えられる。
  * **エラー**: §3.5 の共通の検査。割線・傾き・係数の桁あふれは `NonFinite { index }`。
  * **SciPy との差**: 共通の検査による差（重複・未ソートの $x$、非有限の値）のほかに 1 つある。値が複数列のとき、SciPy は閾値に使う最大値を**全部の列をまとめて** 1 つ取る（`xp.max(f12)`）。そのため、ある列の結果がほかの列の値の大きさに左右され、1 列ずつ補間した結果と一致しないことがある。`interpx` は最大値を列ごとに取る。結果が変わるのは、ある列の $f_1 + f_2$ が、その列の最大値の $10^{-9}$ 倍は超えるが、全列の最大値の $10^{-9}$ 倍は超えない点だけである。スカラー値では SciPy と同じである。
* **`CubicHermiteSpline<T: Scalar, V: Value<T> = T>`（M1-11、SciPy の `CubicHermiteSpline`）**: 各データ点の値と 1 階導関数（傾き）を両方与える $C^1$ の 3 次補間。傾きを決める規則を持たず、与えられた傾きをそのまま使う。導関数が分かっているデータ（微分方程式の解など）や、傾きを自分の規則で決めたいときに使う。
  * **構築**: `CubicHermiteSpline::new(x: &[T], y: &[V], dydx: &[V])`。傾きは値と同じ型 `V` で、値が `[T; M]` なら列ごとの傾きである。点は 2 個以上（SciPy と同じ。2 点なら 3 次式 1 本）。
  * **保持するもの**: 3 次の `PPoly`。`as_ppoly()` で貸し出す。評価・`eval_hinted`・バッチは `PPoly` のものである。係数は共通の手順（3 次エルミート形）で作る。傾きを係数の配列の先頭に写して渡すので、構築の確保は、座標の複製、係数、区間の幅の逆数の 3 個である。
  * **与えた値と傾きの再現**: 区間 $i$ の定数項と 1 次の係数は、$y_i$ と $s_i$ そのもの（ビット一致）である。したがって各区間の左端で、値は $y_i$ に、1 階導関数は $s_i$ に等しい。右端（次の区間の左端でもある）では、左の区間の 3 次式の値と導関数が、丸めの範囲で一致する。
  * **値が `[T; M]` のとき**: 各列を独立に補間する。1 列ずつ作った結果とビット同一になる。
  * **領域外**: 既定は `Extend`（§3.3、SciPy と同じ）。`with_extrapolate` で 4 モードのどれにも変えられる。
  * **周期性は検査しない**: `Periodic` は外挿モードであり、この型に境界条件はない。両端で値や傾きが一致しないデータに `Periodic` を付けると、周期の継ぎ目で値や傾きが跳ぶ（SciPy の `extrapolate='periodic'`、`PchipInterpolator`、`Akima1D` と同じ）。`NotPeriodic` を返すのは、周期を境界条件として選んだ `CubicSpline` だけである。
  * **エラー**: §3.5 の検査を、傾きを含めた順で行う。傾きの長さが違えば `LengthMismatch`、傾きに NaN・無限大があれば `NonFinite { index }`（`index` は `dydx` の中の位置）。`NonFinite` はどの配列かを持たないので、座標・値・傾きのどれであるかは検査の順で決まる（座標が全部有限なら値、値も全部有限なら傾き）。配列を区別する variant は足さない。割線・係数の桁あふれは `NonFinite { index }`。
  * **SciPy との差**: 共通の検査による差（重複・未ソートの $x$）だけである。非有限の値と傾きは SciPy も拒否する。
* **`BSpline1D<T: Scalar, V: Value<T> = T>`（M1-12、SciPy の `make_interp_spline`）**: データ点を通る $k$ 次の B-spline。$k - 1$ 階までの導関数が連続である。
  * **構築**: `BSpline1D::new(x: &[T], y: &[V], degree: usize, bc: BoundaryCondition<V>)`。
  * **保持するもの**: `BSpline`（ノットと制御点）。`as_bspline()` で貸し出す。評価・`eval_hinted`・バッチ・導関数・積分は `BSpline` のものである。
  * **次数**: 1 から 15（`BSPLINE_MAX_DEGREE`）。16 以上は `InvalidParameter { name: "degree", reason: TooLarge { max: 15 } }`、0 は `InvalidParameter { name: "degree", reason: TooSmall { min: 1 } }`。
    * **0 次を持たない理由**: SciPy の 0 次（ノット $x_0, \dots, x_{n-1}, x_{n-1}$、制御点 $y_i$）は、最後のデータ点の値を幅ゼロの区間に置く。`BSpline` は幅ゼロの区間を評価しない（上の「定義域と区間の帰属」）ので、右端で $y_{n-1}$ ではなく $y_{n-2}$ を返すことになり、データを通らない。左の点の値を返す階段は `Nearest1D` の `Previous` が持っている（右端では $y_{n-1}$）。
  * **境界条件と次数の組み合わせ**（SciPy の `bc_type`。`"natural"`・`"clamped"` が使えるのは SciPy でも 3 次だけである）:

    | `bc` | 使える次数 | SciPy の指定 |
    | :--- | :--- | :--- |
    | `NotAKnot`（`Ends(NotAKnot, NotAKnot)` も同じ） | 1〜15 | `None` / `"not-a-knot"` |
    | `Periodic` | 1〜15 | `"periodic"` |
    | `Ends` で、どちらかの端が `FirstDerivative` / `SecondDerivative` | 3 だけ | `"natural"`、`"clamped"`、`([(1, a)], [(2, b)])` |

    3 次以外で端の導関数を指定すると `InvalidParameter { name: "bc", reason: Unsupported }`。端ごとに複数の導関数を与える指定（SciPy の `([(1, a), (2, b)], ...)`。5 次なら端ごとに 2 個要る）と、ノット列の指定（SciPy の `t`）は v0.1.0 の対象外（ロードマップ A-10）。
  * **ノット列**（SciPy と同じ。定義域は $[x_0, x_{n-1}]$ で、両端のノットは多重度 $k + 1$）:
    * **1 次**: 内側のノットは $x_1, \dots, x_{n-2}$。制御点は $y_i$ そのもので、連立方程式を解かない（折れ線）。
    * **`NotAKnot`、奇数次**: 内側のノットは、データ点から両端の $(k+1)/2$ 個ずつを除いたもの（$x_{(k+1)/2}, \dots, x_{n-1-(k+1)/2}$）。
    * **`NotAKnot`、偶数次**: 内側のノットは、隣り合うデータ点の中点 $(x_i + x_{i+1})/2$（$n - 1$ 個）から、両端の $k/2$ 個ずつを除いたもの。中点は $x_i/2 + x_{i+1}/2$ で求める（和が桁あふれしない。あふれない範囲では結果は同じ）。
    * **3 次で端の導関数を指定した端**: その端の隣のデータ点もノットに残す（両端とも導関数なら、内側のノットは $x_1, \dots, x_{n-2}$）。`NotAKnot` の端では、隣のデータ点（$x_1$ または $x_{n-2}$）を除く。
    * **`Periodic`、2 次以上**: 定義域の中のノットは、奇数次では $x_0, \dots, x_{n-1}$。偶数次では両端はそのままで、内側の $x_i$（$1 \le i \le n-2$）を $(x_{i-1} + x_i)/2$ に置き換える。その外側へ、両側に $k$ 個ずつ、ノットの間隔を周期的に繰り返して延ばす（SciPy の `_periodic_knots`）。両端のノットは重ねない。
    * **偶数次で、隣り合う 2 点が近すぎる場合**: 偶数次のノットは 2 点の間に置く。2 点が浮動小数点で隣り合うほど近いと、中点がどちらかの点に丸められ、ノットがデータ点に重なる（系が特異になりうる。重なったノットの区間は幅がゼロで、基底関数を評価できない）。ノットが 2 点の真の間に入らない区間があれば、`DuplicatePoints { first: i, second: i + 1 }` を返す（この方式が区別できない 2 点）。`NotAKnot` では、中点がノットに使われない区間も含めて全部の区間を調べる（規則を 1 つにする）。`Periodic` では、ノットを置く区間（最後の区間以外）を調べる。奇数次と 1 次はノットがデータ点そのものなので、起きない。
    * **ノットの桁あふれ**: 周期のノットは定義域の外へ延びるので、データが有限でもノットがあふれることがある。隣り合う 2 点の差があふれる場合も含め、有限でないノットは `NonFinite { index }`（`index` は近いほうのデータ点）。
  * **連立方程式**（2 次以上）: 未知数は制御点。行は、各データ点での補間条件 $\sum_j B_j(x_i)\,c_j = y_i$ と、端の導関数の条件 $\sum_j B_j^{(\nu)}(x_0)\,c_j = v$ である。左端の導関数の行を最初に、右端のものを最後に置く。こう並べると行列は下帯・上帯とも $k$ 以内の帯行列で、上の帯行列 LU で解く。
    * **帯の幅**: 導関数の行がなければ、帯は上下とも $k - 1$ で足りる。1 点は $k + 1$ 個の基底関数の区間にあるが、いちばん外の対角に届くのは最初の行と最後の行だけで、そこの要素は、区間の端で 0 になる基底関数である（漸化式で厳密に 0 になる）。その 2 つを除いて、帯を $k - 1$ と宣言する。導関数の行があるとき（3 次）は $k$ のままにする。消去の手間は帯の幅の 2 乗に比例するので、狭く宣言するほうが速い（`benches/bspline1d.rs`）。
  * **`Periodic` の連立方程式**: 周期的なノット列では、制御点も周期的になる（$c_{j + n - 1} = c_j$）。未知数を $n - 1$ 個に減らし、$x_0, \dots, x_{n-2}$ での補間条件を並べると、巡回帯行列（下帯・上帯とも $\lfloor k/2 \rfloor$）になる。奇数次では $x_i$ がノットなので、区間の $k + 1$ 個の基底関数のうち最後の 1 個が厳密に 0 で、これを除くと帯が対称になる（3 次なら巡回三重対角）。$n - 1 > 2\lfloor k/2 \rfloor$ なら上の Woodbury の公式で、それ以下なら帯が一周して重なるので、重なる要素を足した密行列として解く。
    * SciPy は $n \le k$ のとき、導関数の連続の条件を並べた別の密行列（$n + k - 1$ 元）を解く。解は同じである（M1-12 の PR で、1〜15 次、2 点以上の全部の組み合わせを SciPy と照合する）。解き方を 1 つにするため、ここでは上の形だけを使う。
  * **点数の下限**: 補間条件と導関数の条件を合わせて、$k$ 次多項式 1 本を決める $k + 1$ 個に足りること。すなわち $n \ge k + 1 - (\text{導関数を指定した端の数})$、かつ 2 点以上。`NotAKnot` は $k + 1$ 点（ちょうどなら $k$ 次多項式 1 本）、3 次で両端とも導関数なら 2 点、`Periodic` は次数によらず 2 点（2 点なら定数）。足りなければ `InsufficientPoints { required, actual }`。
    * `CubicSpline` は 3 点の not-a-knot を放物線、2 点を直線にする（SciPy の `CubicSpline` と同じ）。`BSpline1D` は SciPy の `make_interp_spline` と同じく、これらを受け付けない。
  * **3 次の `BSpline1D` と `CubicSpline`**: 同じデータと同じ境界条件（上の点数を満たすもの）で、同じ曲線になる（丸めの範囲で）。違うのは表現（B-spline と区分多項式）と解き方（帯行列 LU と Thomas 法）で、`CubicSpline` のほうが構築が速い。片端だけ `NotAKnot` の組み合わせは、SciPy の `make_interp_spline` にはないが、`CubicSpline` と同じ `BoundaryCondition` を同じ意味で受け付けるために持つ。
  * **周期性の検査**: `Periodic` では、`CubicSpline` と同じ許容誤差で $y_0$ と $y_{n-1}$ の一致を確かめ、合わなければ `NotPeriodic`。
  * **値が `[T; M]` のとき**: 行列は列によらないので、各列を独立に解いたものとビット同一になる。
  * **領域外**: 既定は `Extend`、境界条件が `Periodic` なら `Periodic`（§3.3）。`with_extrapolate` で 4 モードのどれにも変えられる。
  * **エラー**: 検査の順は、§3.5 の段 1（長さ）→ 次数（`InvalidParameter`）→ 境界条件と次数の組み合わせ（`InvalidParameter`）→ 段 2（点数。下限は上のとおり）→ 段 3・4（非有限、座標の順）→ 境界条件の導関数の値（`InvalidParameter { name: "bc", reason: NonFinite }`）→ 周期性（`NotPeriodic`）→ ノット（偶数次の `DuplicatePoints`、桁あふれの `NonFinite`）→ 求解。ソルバーの破綻と、制御点の桁あふれは `NonFinite { index }`（`index` は該当する行。`CubicSpline` と同じ扱い）。データが有限でも、座標どうしの差が桁あふれするほど離れていると基底関数を評価できず、ソルバーの破綻としてこのエラーになる。
    * **行列は正則か**: `NotAKnot` と端の導関数の行列は、ノットとデータ点が Schoenberg–Whitney の条件を満たすので正則である。Woodbury の公式で解く帯行列 $A$ も、選点行列の正方の部分行列で、同じ条件を満たす。周期の行列全体（$I + V^{\mathsf T} Z$ が正則であること）は、SciPy も前提にしているが、ここでは証明を持たない。2〜8 次、2〜13 点、隣り合う間隔の比が $10^3$ までのランダムなデータ 2 万通りで条件数を調べ、特異に近づく例は見つからなかった（最悪で 8 次の $2.5 \times 10^9$。間隔の比とともに滑らかに大きくなる）。特異な例が見つかったら、`SolverFailed(SolverError::Singular)`（§4。M3 で入れる予定の variant）をここで返すように改める。
  * **SciPy との差**: 共通の検査による差（重複・未ソートの $x$、非有限の値）のほかに、次がある。

    | 項目 | SciPy 1.18 | `interpx` |
    | :--- | :--- | :--- |
    | 0 次 | 階段を返す | `InvalidParameter`（上の理由） |
    | 1 次で `Periodic` | 外挿は `True`（周期にならない） | 外挿は `Periodic`（2 次以上と同じ） |
    | 1 次で 1 点 | 受け付ける（定義域の幅が 0） | `InsufficientPoints`（2 点以上） |
    | 3 次で片端だけ not-a-knot | `make_interp_spline` は受け付けない | 受け付ける（`CubicSpline` と同じ曲線） |
    | 偶数次で、隣り合う 2 点の中点がどちらかの点に丸められる | 特異な行列の例外、または NaN | `DuplicatePoints` |
    | 偶数次の中点の和が桁あふれする | ノットが無限大になる | あふれない形で求める |
    | 端ごとに複数の導関数、ノット列の指定 | 受け付ける | 持たない（A-10） |
* **`BSplineCurve<T: Scalar, const D: usize>`（M1-13）**: $D$ 次元の点列を通るパラメトリック曲線 $u \mapsto [T; D]$。値の型が `[T; D]` の `BSpline1D` に、パラメータの決め方を足したものである。
  * **`BSpline1D` との関係**: SciPy に曲線のクラスはなく、`make_interp_spline(t, points)` に 2 次元の `y` を渡す。これは `BSpline1D<T, [T; D]>` と同じで、各点のパラメータ $t_i$ を呼び出し側が用意する。曲線では点列だけがあって $t_i$ がないのが普通なので、`BSplineCurve` は点列から $t_i$ を決める。内部に `BSpline1D<T, [T; D]>` を持ち、次数・境界条件・ノット列・連立方程式・エラーは `BSpline1D` のものである。
  * **構築**:
    * `BSplineCurve::new(points: &[[T; D]], degree: usize, bc: BoundaryCondition<[T; D]>, parametrization: Parametrization)`: パラメータを点列から決める。
    * `BSplineCurve::with_parameters(t: &[T], points: &[[T; D]], degree: usize, bc: BoundaryCondition<[T; D]>)`: パラメータを与える（SciPy の `make_interp_spline(t, points)`）。結果は同じ引数の `BSpline1D` とビット同一である。
  * **`Parametrization`**（`#[non_exhaustive]`。既定値は持たせず、引数で明示する）: 隣り合う点の距離を $d_i = \lVert p_{i+1} - p_i \rVert_2$ として、

    | 値 | $t_i$（正規化の前） | 性質 |
    | :--- | :--- | :--- |
    | `Uniform` | $i$ | 点の間隔を見ない。間隔が不揃いな点列では波打つ |
    | `ChordLength` | $\sum_{j < i} d_j$ | 曲線に沿った長さの近似。FITPACK の `splprep` の既定と同じ |
    | `Centripetal` | $\sum_{j < i} \sqrt{d_j}$ | 尖った箇所での行き過ぎが少ない |

    どれも最後の値で割って $[0, 1]$ に正規化する（$t_0 = 0$、$t_{n-1} = 1$）。和は前から順に足す。
  * **距離の計算**: $d_i$ は、差の成分のうち絶対値が最大のもので割ってから 2 乗和を取る（途中の 2 乗が桁あふれ・桁落ちしない）。したがって $d_i = 0$ になるのは 2 点が完全に同じときだけである。
  * **保持するもの**: `BSpline1D<T, [T; D]>` と、各データ点のパラメータ。`parameters() -> &[T]` で返す（自動で決めたとき、どの $u$ がどの点に当たるかを知る手段がほかにない）。`as_bspline()` で `BSpline<T, [T; D]>` を貸し出す。接線（1 階導関数）や 2 階導関数は `as_bspline().eval_derivative(u, nu)` で取り、曲線の型に同じ関数を重ねて持たせない。
  * **閉じた曲線**: `BoundaryCondition::Periodic`。最初と最後の点が（`BSpline1D` と同じ許容誤差で）一致していなければ `NotPeriodic`。
  * **領域外**: `BSpline1D` と同じ（既定は `Extend`、`Periodic` のときは `Periodic`）。パラメータの定義域は、`new` では $[0, 1]$、`with_parameters` では $[t_0, t_{n-1}]$。
  * **エラー**: `new` の検査の順は、次数 → 境界条件と次数の組み合わせ → 点数（下限は `BSpline1D` と同じ）→ 点の非有限（`NonFinite { index }`）→ パラメータ → 以降は `BSpline1D`。
    * `ChordLength` と `Centripetal` で隣り合う 2 点が同じなら `DuplicatePoints { first: i, second: i + 1 }`。`Uniform` はパラメータが重ならないので、同じ点が続いても通す。
    * 2 点が同じでなくても、全体の長さに比べて近すぎると、正規化したパラメータが丸めで重なることがある。これも `DuplicatePoints` になる（`BSpline1D` の座標の検査による）。
    * 点が有限でも、差や和が桁あふれするほど離れていれば、パラメータが有限でなくなり `NonFinite { index }`。
    * `with_parameters` のエラーは `BSpline1D::new` と同じ。
  * **持たないもの**: 弧長、弧長による再パラメータ化、最近点の探索。SciPy に対応する機能がなく、この設計の対象外である。平滑化つきの当てはめ（`splprep` の $s > 0$）は FITPACK 族（ロードマップ A-7）。
  * **SciPy との差**: `with_parameters` は `BSpline1D` と同じ（§6.1 の表）。`new` に当たる関数は `make_interp_spline` にはない。パラメータの式は上のとおりで、golden では生成スクリプトが同じ式で $t_i$ を計算して `make_interp_spline` に渡す。

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
* **型**: `RegularGridInterpolator<T: Scalar, const D: usize, V: Value<T> = T>`。$D$ 次元の直交グリッド（軸ごとの座標の直積）の上の値を補間する。`Point = [T; D]`、`Output = V`。
* **構築**: `new(axes: [&[T]; D], values: &[V], method: GridMethod)` と、所有権を受け取る `from_owned(axes: [Vec<T>; D], values: Vec<V>, method: GridMethod)`。
  * `values` は **C 順（row-major、最後の軸が最速）** の平坦な配列（`ndarray` の標準レイアウト・numpy と一致）。点 $(i_0, \dots, i_{D-1})$ の値は `values[((i_0 n_1 + i_1) n_2 + \dots) + i_{D-1}]`。
  * 形は軸の長さ $n_d$ から決まるので、引数に取らない。`values.len()` が $\prod_d n_d$ と違えば `LengthMismatch`。`shape() -> [usize; D]` で返す。
* **`GridMethod`**（`#[non_exhaustive]`）: 構築時に決め、あとから変えない（SciPy は評価のたびに変えられるが、スプラインの方法は構築時に係数を解く）。M2-1 は `Nearest` と `Linear`、M2-2 が `Cubic`・`Quintic`、M2-3 が `Pchip` を足す。SciPy の `slinear` に当たる variant は持たない（下の「`slinear`」）。
* **軸**: 狭義の昇順か、狭義の降順（SciPy と同じ）。軸あたりの点数は、`Nearest` と `Linear` が 2 点以上（グリッドは各次元に幅を持つ）、`Cubic` が 4 点以上、`Quintic` が 6 点以上（次数 + 1。SciPy と同じ）、`Pchip` が 2 点以上（1 次元の `PchipInterpolator` と同じ。SciPy のグリッドは 4 点以上を求める）。
  * 降順の軸は、構築時に昇順へ並べ替えて持つ（軸の座標を逆順にし、値をその軸に沿って逆順にする）。評価はいつも昇順のグリッドで行う。`axes()` と `values()` は並べ替えたあとのものを返す。
* **区間と領域外（軸ごと）**: 各軸で、座標 $x_d$ の区間 $i_d$（$g_{i_d} \le x_d < g_{i_d + 1}$、上端は最後の区間に含める）を二分探索で求める。領域外の扱いは `Extrapolate` の 4 モードで、**全部の軸に同じモードを、軸ごとに独立に**当てる。規則は `PPoly`（§6.1）と同じである（区間探索のコードを共用する）。既定は `Forbid`。

  | モード | 軸ごとの扱い | `Linear` の結果 | `Nearest` の結果 |
  | :--- | :--- | :--- | :--- |
  | `Forbid` | どれかの座標が領域外なら値なし | `None` | `None` |
  | `Extend` | 端の区間を使い、座標はそのまま | 端のセルの式の延長（SciPy の `fill_value=None`） | 端の点の値（SciPy と同じ） |
  | `Clamp` | 座標を $[g_0, g_{n-1}]$ に切り詰める | 箱の表面の値 | `Extend` と同じ |
  | `Periodic` | 座標を $[g_0, g_{n-1})$ に折り返す | 折り返した点の値 | 同じ |

  `Periodic` は値の周期性を検査しない（外挿モードであり、境界条件ではない）。NaN を含む点は `None`。無限大の座標は `PPoly` と同じ（`Forbid`・`Periodic` で `None`、`Clamp` は端、`Extend` は式をそのまま評価する）。軸ごとに別のモードを指定する形は持たない（SciPy にもない）。
* **`Linear`**: 軸 $d$ の正規化した距離を $y_d = (x_d - g_{i_d}) / (g_{i_d + 1} - g_{i_d})$ として、セルの $2^D$ 個の角の値を、重み $\prod_d w_d$（角が下側なら $w_d = 1 - y_d$、上側なら $y_d$）で足す。角は、軸 0 をいちばん外側にして数える順に足す（SciPy と同じ順）。確保しない。
* **`Nearest`**: 軸ごとに、$y_d \le 0.5$ なら下側の点、そうでなければ上側の点を取り、その格子点の値を返す。ちょうど中間は下側である（SciPy と同じ）。降順で渡した軸では、昇順に並べ替えたあとの下側になる（SciPy も並べ替えてから判定する）。
* **ヒント**: `eval_hinted(&[T; D], &mut [Hint; D])`。軸ごとにヒントを持ち、結果は `eval` とビット同一。`eval_batch_into` は、点を 64 点の区切りごとに見て、全部の軸で一方向に、1 点あたり 1 セル以内で進んでいるとき（経路や、画像の行）だけヒントを使う。そうでない並びでは点ごとに探索する（順不同の点ではヒントのほうが遅い。`PPoly` のバッチと同じ判定で、`benches/grid.rs` で測る）。どちらでも結果は `eval` とビット同一である。
* **等間隔の軸**: 区間を割り算で求める近道は入れない（丸めで二分探索と違う区間になりうる）。入れるなら、測ってから別の行で決める。
* **エラー**: 検査の順は §3.5 の段に合わせる。
  1. `LengthMismatch`（`values.len()` と $\prod_d n_d$。積が `usize` に収まらなければ、どんな配列とも一致しないので同じエラー）
  2. `GridAxis { axis, fault: InsufficientPoints { required: 2, actual } }`（軸 0 から順に）
  3. `GridAxis { axis, fault: NonFinite { index } }`（軸 0 から順に）、次に値の `NonFinite { index }`（`index` は平坦な配列での位置）
  4. 軸の向き（軸 0 から順に）: 最初の 2 点が等しければ `DuplicatePoints { first: 0, second: 1 }`。それ以降は、隣と等しければ `DuplicatePoints`、最初の 2 点が決めた向きに逆らえば `NotMonotonic { index }`（逆らった点の位置）。
  * 位置は、渡された順（並べ替える前）のものである。軸の中の位置も、値の平坦な配列での位置も同じで、構築の途中（スプラインの係数、`Pchip` の傾き、軸の選点行列）で見つかった桁あふれの位置も、降順の軸の分を写し戻して返す。
* **SciPy との差**:

  | 項目 | SciPy 1.18 | `interpx` |
  | :--- | :--- | :--- |
  | 領域外の既定 | 例外（`bounds_error=True`） | `None`（`Forbid`）。`eval_or(point, fill)` が `fill_value` に当たる |
  | 値の NaN | 受け付ける（NaN に触れた評価が NaN になる） | `NonFinite` で拒む（§3.5。スプラインの方法は NaN を含む系を解けないので、方法によらず 1 つの規則にする）。欠損のあるグリッドは、欠損を埋めてから渡す |
  | 方法 | 評価のたびに変えられる | 構築時に決める |
  | 1 点だけの軸 | `nearest` と `linear` は受け付ける | 2 点以上 |
  | `slinear` | `linear` と別の名前 | `Linear` に当たる。variant は持たない |
  | `cubic`・`quintic` の係数 | 既定は反復解法 | 直接法。`*_legacy` と一致する |
  | 偏導関数（`nu`） | `cubic`・`quintic` などで取れる | 持たない（A-3） |
  | `Clamp`・`Periodic` | ない | ある |
  | 2 次元の `linear` | 専用の経路（足す順が違う） | どの次元も同じ式。差は丸めの範囲 |
| `pchip` の軸の点数 | 4 点以上 | 2 点以上（1 次元の `PchipInterpolator` と同じ） |
| `pchip` の評価 | 評価のたびに、列全体の 1 次元の補間を作る | 最後の軸の傾きは構築時に求め、残りは点のまわりの値だけから求める。結果は同じ |
* **`slinear`**: SciPy の `slinear` は 1 次のテンソル積スプラインで、`linear` と同じ関数である（領域外の延長も同じ。違うのは計算の経路と丸めだけで、不均一な 2 次元のグリッドで $7 \times 10^{-16}$。M2-2 の設計時の実測）。同じ関数に名前を 2 つ持たせず、`GridMethod::Linear` がそれに当たる。golden の `slinear` のケースは `Linear` と照合する。
* **`Cubic` / `Quintic`（M2-2）**: 各軸について $k$ 次（3 または 5）の、テンソル積 B-spline。SciPy の `make_ndbspl` が作る `NdBSpline` に当たる。
  * **ノット列**: 軸ごとに、`BSpline1D` の not-a-knot と同じ（§6.1。SciPy の `_not_a_knot`）。奇数次だけなので、内側のノットはデータ点から両端の $(k+1)/2$ 個ずつを除いたものである。
  * **係数**: 構築時に 1 度だけ求める。軸ごとの補間演算子は可換なので、軸 0 から順に、その軸に沿った全部の列に 1 次元の補間（選点行列の求解）を当てればよい。選点行列は軸ごとに 1 回だけ分解し（M1-12 の帯行列 LU）、全部の列を同じ分解で解く。`BSpline1D` の「ノット列を作り、選点行列を分解する」部分を crate 内の部品にして共用する。
  * **保持するもの**: 係数（グリッドと同じ形）、軸ごとのノット列、元の値（`values()` が返す）。メモリはグリッド点数の 2 倍に比例する。`Nearest` と `Linear` は値だけを持つ。
  * **評価**: 軸ごとに、区間で 0 でない $k + 1$ 個の基底関数の値をスタック上に求め（`BSpline` と同じ漸化式）、$(k+1)^D$ 個の係数を重みの積で足す。確保しない。区間と領域外は上の表のとおりで、`Extend` は端のセルの多項式の延長である（SciPy の `fill_value=None` と同じ）。
  * **`Periodic`**: 座標を折り返すだけで、スプラインを周期的にはしない（周期の境界条件ではない）。継ぎ目で値は折り返した点のものになり、導関数は一般に跳ぶ。周期の境界条件を持つグリッドのスプラインは持たない（SciPy にもない）。
  * **導関数**: SciPy はこの方法に限り、評価時に偏導関数の階数（`nu`）を取れる。v0.1.0 では持たない（ロードマップ A-3）。
  * **エラー**: M2-1 と同じ順。軸の点数の下限が次数 + 1 になる。ある軸の選点行列を分解できないときは `GridAxis { axis, fault: NonFinite { index } }`（`index` は行。軸の選点行列は正則なので、起きるのは、座標の差が桁あふれするほど離れた軸で要素が有限でなくなる場合だけである。原因が軸にあるので、どの軸かを返す）。係数の桁あふれは `NonFinite { index }`（`index` は係数の平坦な配列での位置）。
* **SciPy 1.18 の既定値との差（`cubic` / `quintic`）**: SciPy はこの 2 メソッドの係数を、既定では反復解法（`gcrotmk`）で求める。その結果はテンソル積スプラインから、滑らかな関数の値では等間隔の 2D グリッドで約 $10^{-5}$、非等間隔や 3D のグリッドでは最大約 $10^{-2}$ ずれ（golden 作成時の実測）、乱数の値では値と同じ大きさまでずれる（不均一な 7×9 のグリッドで、`cubic` が 0.67、`quintic` が 6.5。M2-2 の設計時の実測）。`interpx` は係数を直接法で求めるので、SciPy の既定の出力とは一致しない。
  * **一致する相手**: SciPy の `cubic_legacy` / `quintic_legacy`（軸ごとに 1 次元のスプラインを当てる実装）と、`make_interp_spline` を軸ごとに当てて作った `NdBSpline` である。どちらも同じテンソル積スプラインで、互いに丸めの範囲で一致する。
  * **`spsolve` を基準にしない理由**: `solver=scipy.sparse.linalg.spsolve` も同じスプラインを与えるが、全部の軸をまとめた 1 つの疎行列を解くので、条件数が軸ごとの条件数の積になる。厳密な有理数計算と比べると（不均一な 7×9、乱数の値）、係数の誤差は、最大の係数に対して、軸ごとの解法が $8 \times 10^{-16}$（`cubic`）と $1 \times 10^{-14}$（`quintic`）、`spsolve` が $2 \times 10^{-15}$ と $5 \times 10^{-13}$ だった。`quintic` では `spsolve` のほうが 40 倍ほど粗い。
  * golden の `expected` は `*_legacy` の値とし、既定の値は `expected_scipy_default` として併記する（§8.4.1）。この差は rustdoc にも書く。
* **`Pchip`（M2-3）**: 1 次元の PCHIP（§6.1 の `PchipInterpolator`）を、軸ごとに順に当てる。SciPy の `pchip` に当たる。
  * **軸の順**: 最後の軸から軸 0 へ向かう。まず、最後の軸に沿った列のそれぞれを 1 次元の PCHIP で補間して、座標 $x_{D-1}$ での値を求める。軸が 1 本減ったその値を、軸 $D-2$ に沿って同じように補間し、軸 0 まで続ける。PCHIP の傾きは値に対して線形でないので、軸の順を変えると別の補間になる（乱数の値を置いた不均一な 6×7 のグリッドで、軸を入れ替えると値が $6.5 \times 10^{-3}$ 変わる。設計時の実測）。SciPy と同じ順にする。
  * **点のまわりだけで決まる**: 1 次元の PCHIP は、区間 $[x_i, x_{i+1}]$ の上では $y_{i-1}, \dots, y_{i+2}$ の 4 点だけで決まる（端の区間では 3 点、2 点の軸では 2 点）。両端の傾きが隣の点までしか見ないからである。だから 1 点の評価に要るのは、軸ごとに多くて 4 点、合わせて多くて $4^D$ 個の値である。SciPy は評価のたびに列全体の補間を作るが、同じ結果がこの値だけから求まる。
  * **保持するもの**: 値と、最後の軸に沿った傾き（グリッドと同じ形）。メモリはグリッド点数の 2 倍に比例する（`Cubic` と同じ）。最初の畳み込み（最後の軸に沿った $4^{D-1}$ 本）は元の値に対するものなので、その傾きは評価する点によらず、構築時に求めておける。それより前の軸の畳み込みは、畳み込んだ値が点ごとに違うので、評価時に求める。
    * 決めた根拠（設計時の試作の実測。`benches/grid.rs` と同じグリッドで、1000 点の評価の時間。交互に 3 回ずつ走らせた最良値。経路 / 順不同）:

      | | 2 次元（200 点 × 2） | 3 次元（40 点 × 3） | 4 次元（10 点 × 4） |
      | :--- | :--- | :--- | :--- |
      | 傾きを持たない（値だけ） | 90.8 µs / 93.6 µs | 267 µs / 294 µs | 843 µs / 903 µs |
      | 最後の軸の傾きを持つ | 56.3 µs / 59.3 µs | 145 µs / 154 µs | 427 µs / 457 µs |
      | 比 | 0.62 / 0.63 | 0.54 / 0.52 | 0.51 / 0.51 |

      持たない版は、$4^{D-1}$ 本の列のそれぞれで、評価のたびに傾きを 2 つ求め直す。
    * 最後の軸の 3 次式の係数（区間あたり 4 個）まで持つ形は採らない。メモリが値の 5 倍になる。減らせるのは、最後の軸の 3 次式の係数を作る分だけで、最後の軸の計算は全部でも評価の 1/3 ほどである（3 次元で 156 ns のうち約 53 ns）。
  * **評価**: 軸ごとに区間を求め（上の「区間と領域外」）、$4^{D-1}$ 本の列で、保持した傾きと値から区間の 3 次式を作って評価する。その結果を、軸 $D-2$ から軸 0 へ、軸ごとに 4 個（端では 3 個か 2 個）ずつ畳み込む。畳み込みでは、その 4 個の値から割線・両端の傾き・3 次式を求める。作業用の配列はスタックに置き、確保しない。
    * 傾きの式、区間の 3 次式、その評価は、`PchipInterpolator` のものを crate 内の部品にして共用する。結果は、`PchipInterpolator` を同じ順に軸ごとに当てたものとビット同一である（検査もそう書く）。
    * 傾きの式は SciPy の形のままにする（重みつき調和平均を割り算 4 回で求め、符号の判定で分岐する）。割り算を 2 回に減らした形と、分岐をなくして選択にした形を試したが、どの次元でも速い形はなかった（設計時の実測。順不同の 1000 点、3 回の最良値）:

      | 傾きの式 | 2 次元 | 3 次元 | 4 次元 |
      | :--- | :--- | :--- | :--- |
      | 分岐、割り算 4 回（SciPy の形） | 61.5 µs | 171 µs | 480 µs |
      | 分岐、割り算 2 回 | 56.2 µs | 157 µs | 457 µs |
      | 選択、割り算 4 回 | 52.4 µs | 172 µs | 559 µs |
      | 選択、割り算 2 回 | 51.9 µs | 163 µs | 521 µs |

      選択にすると、傾きが 0 になる点（乱数の値では約半分）でも割り算をするので、列の多い 4 次元で遅くなる。割り算 2 回の形の差（5〜9%）は、同じ版を繰り返したときのばらつき（同じ表の測定で 10% 以上）の内である。SciPy と演算の順まで同じ形を、変える理由がない。
    * 滑らかな値では傾きが 0 にならず、どの点でも割り算をする。その場合も測った（M2-3 の実測。`benches/grid.rs` の `pchip_smooth`、順不同の 1000 点、3 回の最良値）。割り算 4 回の形が 55.7 µs / 163 µs / 522 µs、3 回の形（$(w_1 + w_2) / (w_1/m_{i-1} + w_2/m_i)$）が 55.1 µs / 151 µs / 484 µs、2 回の形が 56.1 µs / 158 µs / 489 µs で、差は 7% 以内である。構築（400×400）でも 7% 以内だった。式はそのままにする。
    * 構築では、区間の幅の逆数を軸ごとに 1 回だけ求める（列ごとに求め直すと、傾きの割り算と同じ数の割り算が増える）。400×400 のグリッドで 2.09 ms が 1.65 ms になった。同じ部品を使う `PchipInterpolator` の構築も 9〜19% 速くなった（1000 点で 9.3 µs が 7.5 µs）。
  * **形を保つ**: 領域内では、結果はその点を含むセルの $2^D$ 個の角の値の、最小と最大の間にある（丸めを除く）。1 次元の PCHIP が区間の両端の値の間に収まり、畳み込みのたびにそれが保たれるからである。値が平らなところでは平らになる。`Cubic` と `Quintic` にはない性質である。領域外へ延ばした先（`Extend`）には及ばない。
  * **領域外**: 既定は `Forbid`。`Extend` は、軸ごとに端の区間の 3 次式を延ばす（SciPy の `fill_value=None` と同じ）。`Clamp` と `Periodic` は、座標を写してから同じ式で求める。
  * **桁あふれ**: 構築時に求めるのは最後の軸の傾きだけなので、検査できるのもそこまでである。最後の軸の割線か傾きが有限でなければ `NonFinite { index }`（`index` は値の平坦な配列での位置）。それより前の軸の割線・傾き・係数は評価時に求めるので、そこであふれたときは、結果が有限でない値になる（IEEE 754 のまま返す。検査を足さない）。値の差が `T` の最大値に近いときにだけ起きる。
  * **エラー**: M2-1 と同じ順。軸の点数の下限は 2。そのあとに、上の `NonFinite`。
  * **持たないもの**: 偏導関数（SciPy も `pchip` では取れない）。複素数の値（`Value` にない。SciPy でも非推奨）。

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
* **ゼロアロケーション検査**: `tests/alloc.rs` がカウンティングアロケーターで、`eval` / `eval_or` / `eval_batch_into` の確保回数が 0 であることを検査する（`just test` に含む）。カウンターはプロセス全体で数えるので、このファイルはテストハーネスなし（`harness = false`、`main` 関数）で動かす。ハーネスは自分のスレッドを持ち、計測中にそこで起きた確保も数えられてしまうためである。並列版 `eval_batch_par` は逐次版とビット同一であることをテストする。
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
* **golden の範囲**: 層 A の golden は、補間族ごとに `compare/goldens/` の `oned` / `cubic_spline` / `bspline1d` / `calculus` / `curve` / `grid` / `grid_pchip` / `scattered` / `rbf` に置く（`cubic_spline` は `CubicSpline` の全ケース、`calculus` は `PPoly` と `BSpline` の導関数・積分で、どちらも 1 ファイルの上限に収めるため `oned` から分けてある。`bspline1d` は `make_interp_spline` の全ケースで、M1-12 で同じ理由から分ける。`grid` は 6 つの方法を同じグリッドで比べるケースで、`grid_pchip` は `pchip` のために選んだケース（軸の入れ替え、段と平らな値、2 点と 3 点の軸など）である。`grid` が上限の 75% にあるので、M2-3 で分ける。SciPy のグリッドが受け付けない 2 点と 3 点の軸のケースは、SciPy の `PchipInterpolator` を最後の軸から順に当てた値を基準にする）。入力のみの `datasets`（古典データ・節点配置・敵対的入力）と、解析関数の値の `functions` も置く。1D のケースは値に加えて、1 階・2 階導関数、原始関数、定積分、区分多項式の係数（または B-spline のノットと制御点）を持つ。`calculus` は、導関数と原始関数そのもの（`PPoly` の係数、`BSpline` のノットと制御点）、3 階以上の導関数、2 階の原始関数、領域外にまたがる定積分、周期をまたぐ定積分を持つ。SciPy と意図して変える挙動（§6.1 の「SciPy との差」）は golden に入れず、解析解と性質テストで確かめる。`oned` には M1-7 で、無限大のクエリ（`numpy.interp`）、`fill_value="extrapolate"` の `nearest` / `nearest-up`、2 列の値（`interp1d`）を加える。M2-1 では `grid` に、`Nearest` でちょうど中間の問い合わせ点、2 列の値、3 次元で一部の軸だけ降順のグリッド、不均一なグリッドと 3 次元での `Extend` を加える。`Clamp` と `Periodic` は SciPy に対応がないので、性質テストで確かめる。M1-13 では `curve` に、3 つのパラメータの決め方それぞれ（パラメータそのものも持つ）、閉じた曲線（周期の 2・3・5 次）、4 次元の点、点数が $k + 1$ ちょうどの点列、間隔が極端に不均一な点列を加える。M1-12 では `bspline1d` に、偶数次（2・4 次）と 6 次以上、点数が $k + 1$ ちょうどのデータ、周期の 1〜5 次と点数が次数以下の周期、3 次の端の導関数が 0 でないもの、3 次で片端だけ not-a-knot のもの（基準は SciPy の `CubicSpline`）、2 列の値、間隔が不均一なデータを加える。帯行列ソルバーの golden は `banded` に置く（LAPACK の `gbsv` の解で、`scipy.linalg.lapack.dgbsv` から得る。`make_interp_spline` が呼ぶものである。`scipy.linalg.solve_banded` は、帯が上下 1 本ずつのとき三重対角専用の `gtsv` を呼ぶので、基準にしない。LAPACK が行を入れ替えた列の数も持ち、ピボットの選び方まで照合する。巡回帯行列は密行列の `scipy.linalg.solve` の解。ランダムな帯行列、行の入れ替えが実際に起きる行列、条件の悪い行列、$n = 1$〜$3$、全幅の帯（密行列）、複数列の右辺）。M1-11 では `CubicHermiteSpline` に、2 点、3 点、間隔が不均一なデータ、傾きが全部 0 のデータ、直線上のデータにその直線の傾きを与えたもの、データと合わない傾き（行き過ぎる例）、2 列の値、周期の外挿（両端が合うデータと、合わないデータ）を加える。M1-10 では `Akima1D` の両方式に、2 点、3 点、4 点、直線上のデータと平らな区間（傾きが式で定まらない点）、2 列の値を加える。M1-9 では `PchipInterpolator` の、2 点、3 点、平らな区間と符号の変化を含むデータ、単調なデータ、2 列の値を加える。M1-8 では `CubicSpline` の、2 点（境界条件ごと）、3 点（not-a-knot・周期・混在）、4 点の周期、not-a-knot と導関数の混在、2 列の値、間隔が極端に不均一なデータを加える。三重対角ソルバーの golden は `tridiagonal` に置く（SciPy の `solve_banded` の解。優対角のランダム行列、古典データから作った 4 種の境界条件の行列、巡回行列、条件の悪い行列、$n = 1$〜$3$、複数列の右辺）。グリッドの `cubic` / `quintic` の基準は §6.4 のとおり `cubic_legacy` / `quintic_legacy` の値とする。M2-2 では `grid` に、係数そのもの（`make_interp_spline` を軸ごとに当てたもの）、点数が下限ちょうどの軸、降順の軸、2 列の値を加える。
* **サイズの上限**: golden は 1 ファイル 1 MB、ディレクトリ合計 8 MB を上限とする（どちらも圧縮しない大きさ）。golden は公開パッケージに含める（取ってきたクレートで `cargo test` が通るようにするため）。crates.io の上限 10 MB は圧縮後の大きさで、golden は約 1/5 に圧縮されるので、8 MB は圧縮後で約 1.6 MB である。1 ファイルの上限は差分とレビューのためのもので、超える系統はファイルを分ける。合計が上限の 90% を超えることになる行は、ケースを落とさず、先に上限を決め直す。経緯と、採らなかった案は ADR 0010。

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
