# `interpx` Architecture & Design Specification

> **SciPy-equivalent, all-in-one interpolation library in Pure Rust.**

---

## 1. ビジョンと設計原則 (Vision & Design Principles)

### 1.1 背景と目的
Rust の科学計算エコシステムにおける最大の課題は「断片化（fragmentation）」である。1次元補間、グリッド補間、散乱点補間が個別のクレートに分かれ、API や型システムが統一されていない。特に不規則散乱データの多次元補間（`LinearNDInterpolator`, `griddata`）は実質的な空白地帯となっている。

`interpx` は、**SciPy (`scipy.interpolate`) の補間アルゴリズム要約表を完全網羅するオールインワンのデファクトスタンダード**を目指す。

### 1.2 コア設計原則
1. **Pure Rust (C/Fortran 依存ゼロ)**: `cargo build` のみでどこでもビルド可能（クロスコンパイル対応）。
2. **ゼロコスト抽象化と静的型安全性**: `const D: usize` によるコンパイル時次元最適化とスタック割り当てを優先。
3. **Rust イディオムと高スループットの両立**: `Option<T>` による安全な境界外判定と、メモリ再利用（ゼロアロケーション）＋ `rayon` ネイティブ並列化を提供。
4. **モジュラーかつ包括的 (Batteries-Included but Modular)**: デフォルトで全機能を提供しつつ、オプトアウトにより最小構成での組み込みも可能。

---

## 2. 外部依存と技術スタック (Technology Stack)

| 役割                 | 採用クレート     | 選定理由                                                                            |
| :------------------- | :--------------- | :---------------------------------------------------------------------------------- |
| **幾何計算エンジン** | **`convx`**      | Pure Rust による $N$ 次元 Delaunay 単体分割・隣接トポロジー・凸包判定。             |
| **線形代数エンジン** | **`faer`**       | Pure Rust 最高速の SIMD 線形代数。スプラインの三重対角解法および RBF の密行列分解。 |
| **空間近傍探索**     | **`kiddo`**      | Pure Rust の高速 k-d tree。静的次元 `[T; D]` と SIMD 最適化に対応。                 |
| **スカラー抽象化**   | **`num-traits`** | `Float` による `f32` / `f64` 完全ジェネリクス。                                     |
| **並列処理 (opt)**   | **`rayon`**      | バッチ補間のゼロコスト並列化。                                                      |
| **配列連携 (opt)**   | **`ndarray`**    | 科学計算エコシステムとのシームレスなデータ受け渡し。                                |
| **エラー定義**       | **`thiserror`**  | 包括的な単一エラー型の定義。                                                        |

---

## 3. 型システムとインターフェース (Type System & APIs)

### 3.1 スカラー型と次元数
* スカラー型は `T: num_traits::Float + ...` によるジェネリクス（`f32` / `f64` 両対応）。
* 多次元座標は原則 `[T; D]`（`const D: usize`）で表現し、動的次元（`&[T]`）を許容するラッパーを提供する。

### 3.2 共通トレイト `Interpolator`
```rust
pub trait Interpolator {
    /// 入力クエリ点の型（1D: T, ND: [T; D] 等）
    type Point;
    /// 出力補間値の型（通常は T）
    type Output;

    /// 単一点の補間。範囲外/外挿不可時は None を返す（型安全）
    fn eval(&self, point: &Self::Point) -> Option<Self::Output>;

    /// 範囲外/外挿不可時に指定の固定値を返す
    fn eval_or(&self, point: &Self::Point, fill_value: Self::Output) -> Self::Output {
        self.eval(point).unwrap_or(fill_value)
    }

    /// バッチ評価（Vec を新しく確保して返却）
    fn eval_batch(&self, points: &[Self::Point]) -> Vec<Option<Self::Output>> {
        points.iter().map(|p| self.eval(p)).collect()
    }

    /// ゼロアロケーション・バッチ評価（既存バッファへの書き込み）
    fn eval_batch_into(
        &self,
        points: &[Self::Point],
        out: &mut [Self::Output],
        fill_value: Self::Output,
    ) {
        assert_eq!(points.len(), out.len());
        for (p, dest) in points.iter().zip(out.iter_mut()) {
            *dest = self.eval_or(p, fill_value);
        }
    }
}
```

### 3.3 構造体の所有権とライフタイム
* すべての補間器インスタンスはライフタイムを持たない（**`'static`**）。
* 前処理済みのインデックス、係数配列、トポロジー構造を自前で所有（Owned）する。
* コンストラクタは `::new(&[...])`（借用スライスから内部バッファを構築）と、ゼロコピー移動が可能な `::from_owned(...)` の両方を提供する。

---

## 4. エラーハンドリング仕様 (Error Handling)

包括的な単一のエラー型 `InterpolateError` を採用する。

```rust
#[derive(Debug, thiserror::Error)]
pub enum InterpolateError {
    #[error("Insufficient points: required at least {required}, but got {actual}")]
    InsufficientPoints { required: usize, actual: usize },

    #[error("Dimension mismatch: expected dimension {expected}, got {actual}")]
    DimensionMismatch { expected: usize, actual: usize },

    #[error("Length mismatch: coordinates length {coords_len} does not match values length {values_len}")]
    LengthMismatch { coords_len: usize, values_len: usize },

    #[error("Input coordinates are not strictly monotonic or contain duplicate points")]
    NonMonotonic,

    #[error("Geometric degeneracy in triangulation: {0}")]
    DegenerateGeometry(String),

    #[error("Linear algebra solver failed: {0}")]
    SolverFailed(String),

    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),
}
```

---

## 5. モジュール構造と SciPy カバレッジ (Coverage Matrix)

`scipy.interpolate` の全要約表に対する実装マッピング。

```text
interpx/
├── 1d/               # 1次元補間
│   ├── linear        # numpy.interp / Linear1D
│   ├── cubic_spline  # CubicSpline (2nd derivative continuous)
│   ├── pchip         # PchipInterpolator (Monotone cubic spline)
│   ├── bspline       # make_interp_spline (k-th degree B-spline)
│   └── nearest       # Nearest1D (nearest, previous, next)
├── curve/            # N次元パラメトリック曲線
│   └── nd_curve      # make_interp_spline (N-D y array / knot-based)
├── grid/             # N次元直交グリッド補間 (Regular / Rectilinear)
│   └── regular_grid  # RegularGridInterpolator (linear, nearest, cubic, quintic, pchip)
└── scattered/        # N次元不規則散乱データ補間
    ├── nearest_nd    # NearestNDInterpolator (kiddo k-d tree)
    ├── linear_nd     # LinearNDInterpolator (convx Delaunay + Jump-and-Walk)
    ├── clough_tocher # CloughTocher2DInterpolator (convx 2D macro-element)
    └── rbf           # RBFInterpolator (faer + kiddo + polynomial tails)
```

### SciPy 全機能網羅対応表

| 分類              | アルゴリズム / 補間ルーチン        | `interpx` 実装構造体                        | 内部コアエンジン                         |
| :---------------- | :--------------------------------- | :------------------------------------------ | :--------------------------------------- |
| **1D**            | Linear (`numpy.interp`)            | `interpx::1d::Linear1D`                     | `PPoly` (ホーナー法)                     |
| **1D**            | Cubic Spline (`CubicSpline`)       | `interpx::1d::CubicSpline`                  | `faer` (三重対角) + `PPoly`              |
| **1D**            | Monotone Cubic Spline (`Pchip`)    | `interpx::1d::PchipInterpolator`            | `PPoly` (単調性制限導関数)               |
| **1D**            | Non-cubic Spline (任意次数 $k$)    | `interpx::1d::BSpline1D`                    | `BSpline` (de Boor)                      |
| **1D**            | Nearest (`interp1d`)               | `interpx::1d::Nearest1D`                    | 二分探索 (Hinted search)                 |
| **N-D Curve**     | Spline N-D Curve                   | `interpx::curve::BSplineCurve<T, D>`        | `BSpline` (ベクトル制御点)               |
| **N-D Grid**      | Nearest Grid                       | `interpx::grid::RegularGridInterpolator`    | インデックス四捨五入                     |
| **N-D Grid**      | Linear Grid (Multilinear)          | `interpx::grid::RegularGridInterpolator`    | テンソル積・次元分離線形結合             |
| **N-D Grid**      | Cubic / Quintic Spline Grid        | `interpx::grid::RegularGridInterpolator`    | 各軸 1D Spline 畳み込み                  |
| **N-D Grid**      | Monotone (PCHIP) Grid              | `interpx::grid::RegularGridInterpolator`    | 各軸 PCHIP 畳み込み                      |
| **N-D Scattered** | Nearest ND (`griddata`)            | `interpx::scattered::NearestNDInterpolator` | `kiddo` (k-d tree)                       |
| **N-D Scattered** | Linear ND (`LinearNDInterpolator`) | `interpx::scattered::LinearNDInterpolator`  | `convx` (Delaunay) + Jump-and-Walk       |
| **N-D Scattered** | Cubic 2D (`CloughTocher2D`)        | `interpx::scattered::CloughTocher2D`        | `convx` (2D 分割) + Bernstein-Bézier     |
| **N-D Scattered** | RBF (`RBFInterpolator`)            | `interpx::scattered::RBFInterpolator`       | `faer` (密行列ソルバー) + `kiddo` (局所) |

---

## 6. コアエンジン内部アルゴリズム詳細 (Internal Algorithms)

### 6.1 スプライン基盤（`PPoly` vs `BSpline`）
* **`PPoly` (Piecewise Polynomial)**:
  * 係数行列 $C \in \mathbb{R}^{(k+1) \times (N-1)}$ を保持。
  * `eval(x)` は二分探索で区間 $i$ を特定した後、ホーナー法 $\sum c_j (x - x_i)^j$ により数回の積和算で高速評価。
* **`BSpline`**:
  * ノットベクトル $t$ と制御点 $c$ を保持。de Boor のアルゴリズムにより任意次数 $k$ を安定評価。

### 6.2 散乱データ幾何探索（`convx` 連携）
* **Delaunay 分割**: `convx` にフラットスライス `&[f64]` を渡し、単体頂点配列と隣接配列（`neighbors`）を取得。
* **Point Location (Jump-and-Walk)**:
  1. **Jump**: `kiddo` の k-d tree でクエリ点 $p$ の最近傍入力頂点を特定し、その頂点を含む単体を初期シードとする。
  2. **Walk**: `convx` の単体隣接グラフを辿り、クエリ点方向の面（Facet）を横断しながら所属単体に到達。
  3. **Barycentric Coordinates**: 到達した単体の頂点座標を行列に組み、重心座標ベクトル $w$ を解く。
  4. 凸包の外側（隣接単体が番兵値 `u32::MAX` で外向き）に突き抜けた場合は即座に `None` を返却。

### 6.3 `RBFInterpolator` の数式モデル
* 連立一次方程式：
  $$
  \begin{bmatrix}
  A & P \\
  P^T & 0
  \end{bmatrix}
  \begin{bmatrix}
  w \\
  c
  \end{bmatrix}
  =
  \begin{bmatrix}
  y \\
  0
  \end{bmatrix}
  $$
  * $A_{ij} = \phi(\|x_i - x_j\|) - \lambda \delta_{ij}$（$\phi$ は基底関数、$\lambda$ は平滑化係数）
  * $P$ は指定次数多項式の Vandermonde 行列（定数項・線形項等のドリフト除去）
* `faer` のコレスキー / LDLT / LU 分解で解を前計算。
* `neighbors` 指定時はクエリごとに `kiddo` で近傍 $k$ 点を抽出して局所部分問題を高速に解く。

---

## 7. Cargo Features 構成 (Crate Features)

```toml
[features]
default = ["1d", "curve", "grid", "scattered", "faer", "convx"]

# モジュール分割（オプトアウト用）
1d = []
curve = ["1d"]
grid = ["1d"]
scattered = ["convx", "faer", "kiddo"]

# 外部クレート連携
convx = ["dep:convx"]
faer = ["dep:faer"]
kiddo = ["dep:kiddo"]
ndarray = ["dep:ndarray"]
rayon = ["dep:rayon"]
serde = ["dep:serde"]
```

---

## 8. ロードマップと完了基準 (Release Checklist)

この仕様に基づき、GitHub 上で以下のチェックリストをすべて満たした時点で `v0.1.0` を正式公開する。

- [ ] **Core Architecture**
  - [ ] 共通トレイト `Interpolator` およびバッチ評価メソッドの実装
  - [ ] 包括的エラー型 `InterpolateError` の定義
  - [ ] 共通数学基盤（`PPoly`, `BSpline`）の実装
- [ ] **1D & Curves**
  - [ ] `Linear1D`, `Nearest1D`
  - [ ] `CubicSpline`（自然境界・固定境界・周期境界）
  - [ ] `PchipInterpolator`（単調性保持）
  - [ ] `BSpline1D` (任意次数 $k$)
  - [ ] `BSplineCurve` (N-D curve)
- [ ] **N-D Regular Grid**
  - [ ] `RegularGridInterpolator` (Nearest, Linear)
  - [ ] `RegularGridInterpolator` (Cubic, Quintic, PCHIP スプライン畳み込み)
- [ ] **N-D Scattered**
  - [ ] `NearestNDInterpolator` (`kiddo` 連携)
  - [ ] `LinearNDInterpolator` (`convx` Delaunay + Jump-and-Walk)
  - [ ] `CloughTocher2DInterpolator` (2D 勾配推定 + Bernstein-Bézier パッチ)
  - [ ] `RBFInterpolator` (`faer` ソルバー + 多項式項 + カーネル 7種 + 局所近傍)
- [ ] **Performance & Integrations**
  - [ ] `rayon` 並列バッチ評価ベンチマーク（SciPy との速度比較）
  - [ ] `ndarray` 相互運用機能
  - [ ] SciPy 出力結果との数値一致テストスイート (Integration Tests)

---

この設計仕様書をベースに、リポジトリ作成からスムーズに開発を進めていただけます！