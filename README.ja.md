[English](README.md) | 日本語

# interpx

純粋な Rust で書いた補間ライブラリです。1 変数、曲線、規則格子、散布点を扱います。入力は構築時に検査します。補間器はすべて自分のデータを持ち、スレッドの間で共有できます。

```toml
[dependencies]
interpx = "0.1"
```

フィーチャー:

- `scattered`（既定で有効）: `NearestNDInterpolator`、`LinearNDInterpolator`、`CloughTocher2DInterpolator`、`RBFInterpolator` を使えるようにします。これらだけが使う依存 crate も、このフィーチャーで入ります。
- `rayon`: 一括評価を複数のスレッドで行います（`ParInterpolator`）。結果は直列の一括評価とビット単位で同じです。
- `serde`: すべての補間器と、その構築に使う引数に `Serialize` と `Deserialize` を実装します。

## 例

```rust
use interpx::{BoundaryCondition, CubicSpline, Extrapolate, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.5, 5.0];
    let y = [1.0, 2.7, 5.8, 4.0, 2.2];

    // not-a-knot の 3 次スプライン。
    let spline = CubicSpline::new(&x, &y, BoundaryCondition::NotAKnot)?;
    let at = spline.eval(&2.5);
    assert!(at.is_some());

    // データの外では値なしにし、埋める値は呼び出し時に渡す。
    let spline = spline.with_extrapolate(Extrapolate::Forbid);
    assert_eq!(spline.eval(&6.0), None);
    let mut out = [0.0; 3];
    spline.eval_batch_into(&[0.5, 4.0, 9.0], &mut out, f64::NAN)?;
    assert!(out[2].is_nan());
    Ok(())
}
```

結果を表示する同じプログラム: `cargo run --example one_variable`。ほかの例: `grid`、`scattered`、`parallel`（`--features rayon`）、`save_and_load`（`--features serde`）。

## 使う前に知っておくこと

- **境界条件は引数です。** `CubicSpline::new(x, y, bc)` と `BSpline1D::new(x, y, k, bc)` は境界条件を受け取ります。`BoundaryCondition::default()` は `NotAKnot` です。
- **埋める値は呼び出し時に渡します。** `eval` は `Option` を返します。外挿しない補間器で定義域の外を問い合わせたとき、および NaN を問い合わせたときは `None` です。`eval_or(point, fill)` と `eval_batch_into(points, out, fill)` は、埋める値を受け取ります。
- **定義域の外での振る舞いは、外挿のモード `Extrapolate` で決まります。** `with_extrapolate` で設定します。既定のモードは、補間の方式ごとに決まっています。
  - `Extend`（端の区間の式をそのまま延ばす）: `CubicSpline`、`PchipInterpolator`、`CubicHermiteSpline`、`BSpline1D`、`BSplineCurve`、`PPoly`、`BSpline`。周期境界条件で作ったスプラインは `Periodic` です。
  - `Forbid`（値なし）: `Akima1D`、`Nearest1D`、`RegularGridInterpolator`。
  - `Clamp`（端の値）: `Linear1D`。
- **入力は構築時に検査します。** 1 変数では、座標が渡された順のまま狭義単調増加である必要があります。データはどこでも有限である必要があり、NaN や無限大を含むデータはエラー（`InterpolateError::NonFinite`）になります。データの NaN が補間の値に紛れ込むことはありません。
- **補間の方式は構築時に決めます。** 評価のたびには変えられません。`RBFInterpolator` はカーネルを引数で受け取ります。
- **値は実数です。** 値を固定長の配列 `[T; M]` にすると、`M` 列をまとめて補間できます。
- **補間器はすべて自分のデータを所有し**、`Send + Sync + 'static` です。保持しておくことも、別のスレッドへ渡すことも、共有することもできます。コンストラクタは渡されたものを複製します。データの一部を渡されたまま保持する型（`PPoly`、`BSpline`、`RegularGridInterpolator`、`NearestNDInterpolator`、`LinearNDInterpolator`、`CloughTocher2DInterpolator`、`RBFInterpolator`）には `from_owned` もあり、そのデータを `Vec` で受け取って複製せずに保持します。
- **座標は `f32` か `f64` です**（1 変数、曲線、規則格子）。散布点は `f64` だけです。
- **規則格子の値は C 順です。** 最後の軸が最も速く変わります。
- **NaN と無限大**: データに含まれるとき、構築時のエラーになります。NaN を問い合わせると値はありません（`None`、または埋める値）。無限大の問い合わせは定義域の外として扱います。`Forbid` と `Periodic` では値なし、`Clamp` では端の値です。`Extend` では端の区間の式をその点で評価するので、値が無限大や NaN になることがあります。
- **全点の `RBFInterpolator` は、構築の間 `N` × `N` の行列を持ちます。** 10 838 点で 940 MB です。近傍の点だけを使う設定（`RbfOptions::with_neighbors`）では、代わりに問い合わせごとに小さな系を解きます。`reciprocal_condition` は、全点の系がどれだけ正確に解けたかを示します。
- **評価はメモリを確保しません。** 例外は 3 つあります。1 つ目は `Vec` を返す `eval_batch` と `eval_batch_par` で、`eval_batch_par` はチャンクごとにも `Vec` を 1 つ確保します。2 つ目は `eval_batch_into_par` で、バッチをスレッドに分けるときに `rayon` が確保します。3 つ目は近傍を使う `RBFInterpolator` で、スレッドで最初に評価するときに作業領域を確保し、そのスレッドが持ち続けます。それより近傍の数が多い補間器をそのスレッドで評価すると、作業領域を広げます。
- **0.x のリリース**: マイナーバージョン（0.1 から 0.2）で公開 API が変わることがあります。ここでいう公開 API は、既定のフィーチャーで公開されるものを指します。
- **この crate にまだないもの**: 規則格子と散布点の補間の導関数、区分多項式の根、多変数の区分多項式、重心補間と Krogh 補間、平滑化スプライン、実行時に決める次元と出力の幅、`no_std`、`RBFInterpolator` の系が特異なときに、代わりに疑似逆行列で解く手段。

## 使い方

### `Interpolator`

補間器はすべて `Interpolator` を実装しています。`eval` は 1 点で評価し、`eval_or` は埋める値つきで評価し、`eval_batch` は新しい `Vec` へ、`eval_batch_into` は呼び出し側のバッファへ一括で評価します。自分の型でも、`eval` を実装すれば一括評価のメソッドが使えます。

```rust
use interpx::{InterpolateError, Interpolator};

/// [0, 1] で 10 x、外では値なし。
struct Ramp;

impl Interpolator for Ramp {
    type Point = f64;
    type Output = f64;

    fn eval(&self, x: &f64) -> Option<f64> {
        (0.0..=1.0).contains(x).then(|| 10.0 * x)
    }
}

fn main() -> Result<(), InterpolateError> {
    assert_eq!(Ramp.eval(&0.5), Some(5.0));
    assert_eq!(Ramp.eval_or(&2.0, -1.0), -1.0);
    assert_eq!(Ramp.eval_batch(&[0.5, 2.0]), [Some(5.0), None]);
    let mut out = [0.0; 2];
    Ramp.eval_batch_into(&[0.5, 2.0], &mut out, -1.0)?;
    assert_eq!(out, [5.0, -1.0]);
    Ok(())
}
```

### `Hint`

`eval_hinted` は、直前の問い合わせの区間を `Hint` に覚えておきます。順に並んだ問い合わせは、その区間から 1 歩で自分の区間を見つけます。値は `eval` と同じです。

```rust
use interpx::{Hint, InterpolateError, Interpolator, Linear1D};

fn main() -> Result<(), InterpolateError> {
    let line = Linear1D::new(&[0.0, 1.0, 2.0, 3.0], &[0.0, 10.0, 20.0, 30.0])?;
    let mut hint = Hint::new();
    for x in [0.5, 1.5, 2.5] {
        assert_eq!(line.eval_hinted(x, &mut hint), line.eval(&x));
    }
    Ok(())
}
```

### `Extrapolate`

`Forbid`（値なし）、`Clamp`（近い方の端の値）、`Extend`（端の区間の式を延ばす）、`Periodic`（問い合わせを定義域の中へ折り返す）。

```rust
use interpx::{Extrapolate, InterpolateError, Interpolator, Linear1D};

fn main() -> Result<(), InterpolateError> {
    let line = Linear1D::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 0.0])?;
    assert_eq!(line.extrapolate(), Extrapolate::Clamp);
    assert_eq!(line.eval(&3.0), Some(0.0));
    assert_eq!(line.clone().with_extrapolate(Extrapolate::Forbid).eval(&3.0), None);
    assert_eq!(line.clone().with_extrapolate(Extrapolate::Extend).eval(&3.0), Some(-1.0));
    assert_eq!(line.with_extrapolate(Extrapolate::Periodic).eval(&2.5), Some(0.5));
    Ok(())
}
```

### エラー: `InterpolateError`、`AxisFault`、`DegeneracyKind`、`ParamReason`、`SolverError`

コンストラクタは `InterpolateError` を返し、どの入力のどこが悪いかを示します。`InsufficientPoints`、`LengthMismatch`、`NonFinite`、`NotSorted`、`DuplicatePoints`、`NotPeriodic`、`InvalidParameter`（`ParamReason` つき）、`SolverFailed`（`SolverError` つき）、`AllocationFailed`、`DegenerateGeometry`（散布点を三角形に分けられないとき。`DegeneracyKind` つき）、`GridAxis`（規則格子の 1 本の軸について。`AxisFault` つき）があります。

```rust
use interpx::{
    AxisFault, BSpline1D, BoundaryCondition, DegeneracyKind, GridMethod, InterpolateError,
    Linear1D, LinearNDInterpolator, ParamReason, RBFInterpolator, RbfKernel,
    RegularGridInterpolator, SolverError,
};

fn main() {
    let unsorted = Linear1D::new(&[0.0, 2.0, 1.0], &[0.0, 1.0, 2.0]);
    assert_eq!(unsorted.err(), Some(InterpolateError::NotSorted { index: 2 }));

    let degree = BSpline1D::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 0.0], 16, BoundaryCondition::NotAKnot);
    assert_eq!(
        degree.err(),
        Some(InterpolateError::InvalidParameter { name: "degree", reason: ParamReason::TooLarge { max: 15 } }),
    );

    let axis = RegularGridInterpolator::new([&[0.0, 1.0][..], &[0.0, 0.0]], &[0.0; 4], GridMethod::Linear);
    assert_eq!(
        axis.err(),
        Some(InterpolateError::GridAxis { axis: 1, fault: AxisFault::DuplicatePoints { first: 0, second: 1 } }),
    );

    // 直線上の点と、多項式の項が 2 番目の座標を含むカーネル: その係数が
    // 決まらない。
    let line = [[0.0, 3.0], [1.0, 3.0], [2.0, 3.0], [3.0, 3.0]];
    let singular = RBFInterpolator::new(&line, &[0.0, 1.0, 4.0, 9.0], RbfKernel::ThinPlateSpline);
    assert_eq!(singular.err(), Some(InterpolateError::SolverFailed(SolverError::Singular)));

    // 同じ直線上の点では三角形ができない。
    let flat = LinearNDInterpolator::new(&line, &[0.0, 1.0, 4.0, 9.0]);
    assert_eq!(
        flat.err(),
        Some(InterpolateError::DegenerateGeometry(DegeneracyKind::LowerDimensional { affine_dim: 1 })),
    );
}
```

### `Scalar` と `Value`

`Scalar` は `f32` か `f64` です。`Value<T>` は `T`、または各点に `M` 個の値を持つ `[T; M]` です。`[T; M]` は列ごとに補間され、`M` 個の補間器を別々に作った場合と同じ値になります。

```rust
use interpx::{InterpolateError, Interpolator, Linear1D, Scalar, Value};

fn halfway<T: Scalar, V: Value<T>>(x: &[T], y: &[V]) -> Result<Option<V>, InterpolateError> {
    let line = Linear1D::new(x, y)?;
    Ok(line.eval(&((x[0] + x[x.len() - 1]) / (T::one() + T::one()))))
}

fn main() -> Result<(), InterpolateError> {
    assert_eq!(halfway(&[0.0_f32, 2.0], &[1.0_f32, 3.0])?, Some(2.0));
    assert_eq!(halfway(&[0.0, 2.0], &[[1.0, 10.0], [3.0, 30.0]])?, Some([2.0, 20.0]));
    Ok(())
}
```

### `PPoly`

区分多項式です。区切りの点と、区間ごとの係数（次数の高い順）から作ります。`BSpline` と同じく、`derivative`、`antiderivative`、`eval_derivative`、`integrate` を持ちます。

```rust
use interpx::{InterpolateError, Interpolator, PPoly};

fn main() -> Result<(), InterpolateError> {
    // [0, 1] で x、[1, 2] で 1 + 2(x - 1)。
    let p = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?;
    assert_eq!(p.eval(&1.5), Some(2.0));
    assert_eq!(p.eval_derivative(1.5, 1), Some(2.0));
    assert_eq!(p.integrate(0.0, 2.0), Some(2.5));
    assert_eq!(p.derivative(1)?.eval(&0.5), Some(1.0));
    assert_eq!(p.antiderivative(1)?.eval(&1.0), Some(0.5));
    Ok(())
}
```

### `BSpline` と `BSPLINE_MAX_DEGREE`

ノットと制御点から作る B スプラインです。次数は `BSPLINE_MAX_DEGREE`（15）までです。`PPoly` と相互に変換できます。

```rust
use interpx::{BSpline, InterpolateError, Interpolator, PPoly, BSPLINE_MAX_DEGREE};

fn main() -> Result<(), InterpolateError> {
    // [0, 1] の 2 次式 2 x (1 - x)。
    let b = BSpline::new(&[0.0, 0.0, 0.0, 1.0, 1.0, 1.0], &[0.0, 1.0, 0.0], 2)?;
    assert_eq!(b.eval(&0.5), Some(0.5));
    let area = b.integrate(0.0, 1.0).unwrap_or(f64::NAN);
    assert!((area - 1.0 / 3.0).abs() < 1e-15);

    let p = PPoly::from_bspline(&b)?;
    assert_eq!(p.eval(&0.5), Some(0.5));
    assert!(BSpline::from_ppoly(&p)?.degree() <= BSPLINE_MAX_DEGREE);
    Ok(())
}
```

### `Linear1D`

線形補間です。データ点での値は、データの値と厳密に一致します。

```rust
use interpx::{InterpolateError, Interpolator, Linear1D};

fn main() -> Result<(), InterpolateError> {
    let line = Linear1D::new(&[0.0, 1.0, 3.0], &[0.0, 2.0, 6.0])?;
    assert_eq!(line.eval(&2.0), Some(4.0));
    assert_eq!(line.domain(), (0.0, 3.0));
    Ok(())
}
```

### `Nearest1D` と `NearestKind`

近くのデータ点の値を返します。いちばん近い点（等距離なら左か右）、直前の点、直後の点から選びます。

```rust
use interpx::{InterpolateError, Interpolator, Nearest1D, NearestKind};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0];
    let y = [10.0, 20.0, 30.0];
    assert_eq!(Nearest1D::new(&x, &y, NearestKind::Nearest)?.eval(&0.5), Some(10.0));
    assert_eq!(Nearest1D::new(&x, &y, NearestKind::NearestUp)?.eval(&0.5), Some(20.0));
    assert_eq!(Nearest1D::new(&x, &y, NearestKind::Previous)?.eval(&0.9), Some(10.0));
    assert_eq!(Nearest1D::new(&x, &y, NearestKind::Next)?.eval(&0.1), Some(20.0));
    Ok(())
}
```

### `CubicSpline`、`BoundaryCondition`、`EndCondition`

C² 級の 3 次スプラインです。`BoundaryCondition` は `NotAKnot`、`Periodic`、または両端に `EndCondition`（`NotAKnot`、`FirstDerivative`、`SecondDerivative`）を指定する `Ends` です。`natural()` は両端で 2 階導関数を 0 に、`clamped()` は両端で 1 階導関数を 0 にします。

```rust
use interpx::{BoundaryCondition, CubicSpline, EndCondition, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.0];
    let y = [0.0, 1.0, 0.0, 1.0];
    let natural = CubicSpline::new(&x, &y, BoundaryCondition::natural())?;
    let clamped = CubicSpline::new(&x, &y, BoundaryCondition::clamped())?;
    let slopes = BoundaryCondition::Ends(EndCondition::FirstDerivative(1.0), EndCondition::SecondDerivative(0.0));
    let mixed = CubicSpline::new(&x, &y, slopes)?;
    assert_eq!(natural.eval(&1.0), Some(1.0));
    assert_eq!(clamped.eval(&2.0), Some(0.0));
    assert!(mixed.eval(&1.5).is_some());

    // 周期境界条件: 最初と最後の値が等しい。
    let periodic = CubicSpline::new(&x, &[0.0, 1.0, 0.5, 0.0], BoundaryCondition::Periodic)?;
    assert!(periodic.eval(&4.5).is_some());
    Ok(())
}
```

### `PchipInterpolator`

形を保つ C¹ 級の 3 次補間（PCHIP）です。データが単調な区間では補間も単調で、データの値の範囲を越えません。

```rust
use interpx::{InterpolateError, Interpolator, PchipInterpolator};

fn main() -> Result<(), InterpolateError> {
    let pchip = PchipInterpolator::new(&[0.0, 1.0, 2.0, 3.0], &[0.0, 0.0, 1.0, 1.0])?;
    let v = pchip.eval(&1.5).unwrap_or(f64::NAN);
    assert!((0.0..=1.0).contains(&v));
    Ok(())
}
```

### `Akima1D` と `AkimaMethod`

秋間の C¹ 級の 3 次補間です。重みは、秋間の重みと、その修正版（makima）から選べます。各点の傾きを、近くの 5 点から決めます。

```rust
use interpx::{Akima1D, AkimaMethod, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let makima = Akima1D::new(&x, &y, AkimaMethod::Makima)?;
    // 平らな区間は平らなまま。
    assert_eq!(makima.eval(&0.5), Some(0.0));
    assert!(Akima1D::new(&x, &y, AkimaMethod::Akima)?.eval(&2.5).is_some());
    Ok(())
}
```

### `CubicHermiteSpline`

与えた値と傾きを通る 3 次補間です。

```rust
use interpx::{CubicHermiteSpline, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    // -1、0、1 での値と傾きから x³ を作る。
    let cube = CubicHermiteSpline::new(&[-1.0, 0.0, 1.0], &[-1.0, 0.0, 1.0], &[3.0, 0.0, 3.0])?;
    assert_eq!(cube.eval(&0.5), Some(0.125));
    Ok(())
}
```

### `BSpline1D`

データ点を通る 1〜15 次の B スプラインです。`as_bspline` で、中身の `BSpline` を参照できます。

```rust
use interpx::{BSpline1D, BoundaryCondition, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0];
    let quintic = BSpline1D::new(&x, &y, 5, BoundaryCondition::NotAKnot)?;
    assert_eq!(quintic.degree(), 5);
    let v = quintic.eval(&2.5).unwrap_or(f64::NAN);
    assert!((v - 6.25).abs() < 1e-12);
    assert_eq!(quintic.as_bspline().degree(), 5);
    Ok(())
}
```

### `BSplineCurve` と `Parametrization`

`D` 次元の点列を通るパラメトリック曲線です。`new` は、決め方（`Uniform`、`ChordLength`、`Centripetal`）に従ってパラメータを求めます。`with_parameters` は、パラメータを呼び出し側から受け取ります。

```rust
use interpx::{BSplineCurve, BoundaryCondition, InterpolateError, Interpolator, Parametrization};

fn main() -> Result<(), InterpolateError> {
    let points = [[0.0, 0.0], [1.0, 2.0], [3.0, 3.0], [4.0, 1.0]];
    let curve = BSplineCurve::new(&points, 3, BoundaryCondition::NotAKnot, Parametrization::Centripetal)?;
    let (start, end) = curve.domain();
    assert_eq!(curve.eval(&start), Some(points[0]));
    assert!(curve.eval(&end).is_some());

    let t = [0.0, 1.0, 2.0, 3.0];
    let given = BSplineCurve::with_parameters(&t, &points, 3, BoundaryCondition::NotAKnot)?;
    assert_eq!(given.parameters(), &t);
    Ok(())
}
```

### `RegularGridInterpolator` と `GridMethod`

`D` 次元の直交格子の上の値を補間します。方式は `Nearest`、`Linear`、`Cubic`、`Quintic`、`Pchip` です。軸の間隔は不均一でもよく、降順で与えることもできます。

```rust
use interpx::{Extrapolate, GridMethod, InterpolateError, Interpolator, RegularGridInterpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.0];
    let y = [0.0, 2.0, 4.0, 6.0];
    // f(x, y) = x + y を C 順で: values[i * 4 + j] = x[i] + y[j]。
    let values: Vec<f64> = x.iter().flat_map(|a| y.iter().map(move |b| a + b)).collect();
    let grid = RegularGridInterpolator::new([&x, &y], &values, GridMethod::Linear)?;
    assert_eq!(grid.eval(&[1.5, 3.0]), Some(4.5));
    assert_eq!(grid.shape(), [4, 4]);
    assert_eq!(grid.eval(&[4.0, 0.0]), None);
    assert_eq!(grid.with_extrapolate(Extrapolate::Extend).eval(&[4.0, 0.0]), Some(4.0));

    let cubic = RegularGridInterpolator::new([&x, &y], &values, GridMethod::Cubic)?;
    let v = cubic.eval(&[1.5, 3.0]).unwrap_or(f64::NAN);
    assert!((v - 4.5).abs() < 1e-12);
    Ok(())
}
```

### `NearestNDInterpolator`

`D` 次元で、いちばん近いデータ点の値を返します。`rescaled` は、各座標を同じ範囲にそろえてから距離を測ります。`with_distance_bound` は、ある距離より遠い問い合わせを値なしにします。

```rust
use interpx::{InterpolateError, Interpolator, NearestNDInterpolator};

fn main() -> Result<(), InterpolateError> {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    let nearest = NearestNDInterpolator::new(&points, &[1.0, 2.0, 3.0])?;
    assert_eq!(nearest.eval(&[0.9, 0.2]), Some(2.0));

    let bounded = nearest.with_distance_bound(0.5);
    assert_eq!(bounded.eval(&[5.0, 5.0]), None);

    let rescaled = NearestNDInterpolator::rescaled(&points, &[1.0, 2.0, 3.0])?;
    assert_eq!(rescaled.len(), 3);
    Ok(())
}
```

### `LinearNDInterpolator`

`D` 次元で、点の Delaunay 分割の上で線形に補間します。点の凸包の外では値がありません。`rescaled` は、各座標を同じ範囲にそろえてから三角形分割します。`with_simplex_tolerance` は、単体の内側とみなす範囲を広げます。

```rust
use interpx::{InterpolateError, Interpolator, LinearNDInterpolator};

fn main() -> Result<(), InterpolateError> {
    // f(x, y) = 1 + 2x + 3y at the corners of the unit square.
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
    let linear = LinearNDInterpolator::new(&points, &[1.0, 3.0, 4.0, 6.0])?;
    let v = linear.eval(&[0.25, 0.5]).unwrap_or(f64::NAN);
    assert!((v - 3.0).abs() < 1e-12);
    assert_eq!(linear.eval(&[1.5, 0.5]), None);

    let rescaled = LinearNDInterpolator::rescaled(&points, &[1.0, 3.0, 4.0, 6.0])?;
    assert_eq!(rescaled.simplex_count(), 2);
    Ok(())
}
```

### `CloughTocher2DInterpolator`

平面の点の Delaunay 分割の上で、区分 3 次の補間をします。値も 1 階導関数も連続です。点での傾きは、曲面が辺に沿ってなるべく曲がらないように、データから推定します。推定した傾きは `gradients()` で取り出せます。点の凸包の外では値がありません。`rescaled`、`with_simplex_tolerance`、問い合わせ点を含む三角形の探し方は、`LinearNDInterpolator` と同じです。

```rust
use interpx::{CloughTocher2DInterpolator, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    // 単位正方形のいくつかの点での f(x, y) = x^2 + y。
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.5], [0.3, 0.8]];
    let values: Vec<f64> = points.iter().map(|p| p[0] * p[0] + p[1]).collect();
    let cubic = CloughTocher2DInterpolator::new(&points, &values)?;
    // データ点では、与えた値。
    let v = cubic.eval(&[0.5, 0.5]).unwrap_or(f64::NAN);
    assert!((v - 0.75).abs() < 1e-12);
    assert_eq!(cubic.eval(&[1.5, 0.5]), None);
    assert_eq!(cubic.gradients().len(), points.len());
    Ok(())
}
```

### `RBFInterpolator`、`RbfKernel`、`RbfDegree`、`RbfOptions`

`D` 次元の放射基底関数による補間です。カーネルは 8 つあり、形状パラメータのある 4 つでは `epsilon` を指定します。`RbfOptions` で、平滑化、多項式の項の次数（`RbfDegree`）、近傍の点の数を設定します。

```rust
use interpx::{InterpolateError, Interpolator, RBFInterpolator, RbfDegree, RbfKernel, RbfOptions};

fn main() -> Result<(), InterpolateError> {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.5]];
    let values = [0.0, 1.0, 1.0, 2.0, 1.2];

    let tps = RBFInterpolator::new(&points, &values, RbfKernel::ThinPlateSpline)?;
    let v = tps.eval(&[1.0, 0.0]).unwrap_or(f64::NAN);
    assert!((v - 1.0).abs() < 1e-12);
    assert!(tps.reciprocal_condition().is_some());

    let options = RbfOptions::new()
        .with_smoothing(0.1)
        .with_degree(RbfDegree::Degree(1))
        .with_neighbors(4);
    let gaussian = RbfKernel::Gaussian { epsilon: 2.0 };
    let local = RBFInterpolator::with_options(&points, &values, gaussian, &options)?;
    assert_eq!(local.neighbors(), Some(4));
    assert!(local.eval(&[0.3, 0.6]).is_some());
    Ok(())
}
```

### `ParInterpolator`（フィーチャー `rayon`）

`Sync` な補間器すべてで、`eval_batch_par` と `eval_batch_into_par` が使えます。一括評価は最初の何点かで時間を測り、分ける価値があるときだけ残りを複数のスレッドに分けます。結果は直列の一括評価とビット単位で同じです。

```rust
use interpx::{InterpolateError, Interpolator, Linear1D, ParInterpolator};

fn main() -> Result<(), InterpolateError> {
    let line = Linear1D::new(&[0.0, 1.0], &[0.0, 10.0])?;
    let queries: Vec<f64> = (0..10_000).map(|i| f64::from(i) / 10_000.0).collect();
    assert_eq!(line.eval_batch_par(&queries), line.eval_batch(&queries));
    let mut out = vec![0.0; queries.len()];
    line.eval_batch_into_par(&queries, &mut out, f64::NAN)?;
    Ok(())
}
```

### 保存と読み戻し（フィーチャー `serde`）

補間器はすべて、`serde` のどのフォーマットでも書き出せ、作り直さずに読み戻せます。保存されるのは構築で求めたもの（スプラインならデータではなく係数）と、保存データのバージョンです。読み込みではコンストラクタと同じ検査を行うので、壊れた保存データはエラーになり、パニックする補間器にはなりません。フォーマットがすべての数を元のとおりに戻す場合（バイナリのフォーマット、`float_roundtrip` フィーチャーを有効にした `serde_json`）、読み戻したものは保存したものとビット単位で同じ値を返します。

```rust
use interpx::{BoundaryCondition, CubicSpline, Interpolator};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let spline = CubicSpline::new(&[0.0, 1.0, 2.5, 4.0], &[1.0, 3.0, 2.0, 5.0], BoundaryCondition::NotAKnot)?;
    let text = serde_json::to_string(&spline)?;
    let back: CubicSpline<f64> = serde_json::from_str(&text)?;
    assert_eq!(back.eval(&1.7), spline.eval(&1.7));
    Ok(())
}
```

## 検証

結果は SciPy 1.18 と照らし合わせて検証しています。テストは、SciPy 1.18 が出した値（`Linear1D` だけは NumPy 2.5 の `numpy.interp` の値）と、すべての補間器の値を、決まった許容誤差で比べます。この値はパッケージの `compare/goldens/` に含まれます。SciPy 1.18 が別のものを計算する場合は、結果が次のように違います。

- **規則格子の `Cubic` と `Quintic`**: interpx はスプラインを厳密に解き、SciPy の `method="cubic_legacy"` と `"quintic_legacy"` に丸め誤差の範囲で一致します。SciPy の `cubic` と `quintic` は既定で反復解法を使います。その値とスプラインのずれは、滑らかなデータ（2 次元・3 次元の格子）で値の最大の絶対値の 1e-5〜1e-2 倍、乱数の値では値そのものの大きさに達しました。
- **`CloughTocher2DInterpolator`**: interpx は、残差が右辺の 1e-14 倍を下回るまで傾きを解きます。SciPy は `tol=1e-6` か `maxiter=400` で止まります。比べた 11 組のデータでは、SciPy の値と解のずれは、最大で値の最大の絶対値の 3e-8 倍でした。
- **`LinearNDInterpolator` と `CloughTocher2DInterpolator`**: 同じ円の上に 4 点以上（球の上なら 5 点以上、…）あると、それらが作る区画を、点の番号による規則で分割します。SciPy の三角形分割（Qhull）は別の分け方をします。格子の上の点がこれに当たります。線形補間では、その区画の中で値が違います。Clough-Tocher 補間では、傾きがすべての辺をまとめた 1 つの連立方程式から決まるので、全体で値が違います。座標が同じ点が複数あるとき、interpx は最初に渡された点の値を使います。Qhull はそのうちの 1 つを残しますが、最初の点とは限りません。`LinearNDInterpolator` は 1 次元の点も扱えますが、SciPy は受け付けません。
- **`RBFInterpolator`**: 系は LU 分解で解き、SciPy は対称な不定値の分解を使います。テストで解いた系では、値の差は値の最大の絶対値の 1e-12 倍未満です。SciPy には既定のカーネル（thin plate spline）があり、どのカーネルでも `epsilon` を受け取ります。interpx ではカーネルは引数で、`epsilon` を持つのは形状パラメータのある 4 つのカーネルです。SciPy は、カーネルが求める最小より低い多項式の次数には警告を出し、平滑化があれば重複した点も受け付けます。interpx はどちらも拒否します。有限でない問い合わせは値なしです（SciPy は NaN）。近傍を使うとき、近傍の系が特異になる問い合わせは値なしです（SciPy は例外を出します）。同じ距離に近傍が複数あるとき、interpx は最初に渡された点を選びます。SciPy とは違う点になることがあります。
- **1 変数**: SciPy の `interp1d` は座標を並べ替え、`numpy.interp` は渡された座標を検査せずに使います。interpx は狭義単調増加でない座標を拒否し、2 点以上を求めます（SciPy と NumPy は 1 点でも受け付けます）。`Nearest1D` は NaN の問い合わせに値を返しません（SciPy は最後のデータの値）。`Nearest1D` の `Clamp` と `Extend` は、どの種類でも両側で端の値を返します（SciPy の `"previous"` と `"next"` は片側で NaN を返します）。`Akima1D` は配列の値を列ごとに補間します。SciPy は、重みが消える位置を全部の列をまとめて決めます。
- **`NearestNDInterpolator`**: 有限でない問い合わせは値なしです（SciPy は例外を出します）。座標が同じ点が複数あるとき、interpx は最初に渡された点を選び、SciPy の選び方は木の形によります。`rescaled` は座標ごとの平均を SciPy とは別の順に足すので、距離がちょうど上限にある問い合わせでは判定が違うことがあります。
- **`LinearNDInterpolator`**: 丸め誤差で凸包の面からわずかに外へ出た問い合わせは、すべての単体から探します。重心座標は LU 分解で解きます。SciPy は逆行列を掛けます。条件数が約 1e4 の四面体で 1 次関数を補間すると、誤差は interpx が 1.8e-15、SciPy が 1.8e-12 でした。
- **`BSpline1D`**: 0 次は拒否します（階段関数は `NearestKind::Previous` の `Nearest1D` です）。3 次では、片端に not-a-knot、もう片端に導関数を指定できますが、SciPy は拒否します。1 次の `Periodic` は周期的に外挿しますが、SciPy はしません。偶数次で、2 つのデータ点の中点に置くノットが丸めでデータ点に重なるときは `DuplicatePoints` です（SciPy は特異行列か NaN）。そのまま足すとあふれる 2 点でも、中点はあふれずに求めます（SciPy はノットが無限大になります）。
- **データの NaN と無限大**: 構築時に拒否します。SciPy の `RegularGridInterpolator` は NaN を受け付けてその近くで NaN を返し、`numpy.interp` は無限大を受け付けます。
- **`RegularGridInterpolator` の格子の外**: `Forbid` では値なしです。SciPy は既定で例外を出します（`bounds_error=True`）。
- **`PPoly` と `BSpline`**: `Forbid` のとき、積分の上限か下限が定義域の外にある `integrate` は `None` です（SciPy は、`PPoly` では NaN、`BSpline` では定義域の中の部分の積分）。上限か下限が NaN のときは、どのモードでも `None` です（SciPy は、`PPoly` では例外、`BSpline` では NaN）。`BSpline` で次数より高い階数の導関数は 0 のスプラインで、内側のノットが `degree + 1` 回重なる `BSpline` は区間ごとに微分します。SciPy の `BSpline` は、どちらにも例外を出します。`BSpline::antiderivative` は、ノットによらず定義域の左端で 0 です。SciPy では、最初のノットが `degree + 1` 回重なっているときにそうなります。`BSpline` から作る `PPoly` は、重なったノットで幅 0 の区間を持ちません。SciPy は持ちます。
- **規則格子の `Pchip`** は、軸ごとに 2 点あれば使えます。SciPy は 4 点を求めます。

## 性能

SciPy 1.18 と、ほかの純粋な Rust の補間 crate との比較（時間とメモリ。マシン、バージョン、データ、コマンドつき）: [docs/performance.md](https://github.com/YUKIKEDA/interpx/blob/v0.1.1/docs/performance.md)。

## ライセンス

[Apache License, Version 2.0](LICENSE-APACHE) と [MIT license](LICENSE-MIT) のどちらかを選んで利用できます。
