[English](README.md) | 日本語

# interpx

SciPy 1.18 の `scipy.interpolate` の主なクラスを、純粋な Rust で実装した補間ライブラリです。1 変数、曲線、規則格子、散布点を扱います。SciPy が同じ補間を計算する場合、値は丸め誤差の範囲で SciPy と一致し、テストは SciPy 1.18 が出した値と照合しています。

```toml
[dependencies]
interpx = "0.1"
```

フィーチャー:

- `scattered`（既定で有効）: `NearestNDInterpolator`、`LinearNDInterpolator`、`RBFInterpolator` を使えるようにします。これらだけが使う依存 crate も、このフィーチャーで入ります。
- `rayon`: 一括評価を複数のスレッドで行います（`ParInterpolator`）。結果は直列の一括評価とビット単位で同じです。
- `serde`: すべての補間器と、その構築に使う引数に `Serialize` と `Deserialize` を実装します。

## SciPy から移る方へ

最初の例の前に、SciPy と違う点をまとめます。

- **境界条件は引数です。** `CubicSpline::new(x, y, bc)` と `BSpline1D::new(x, y, k, bc)` は境界条件を受け取ります。SciPy の既定に当たるのは `BoundaryCondition::NotAKnot`（`BoundaryCondition::default()` も同じ）です。
- **構築時の `fill_value` はありません。** `eval` は `Option` を返します。外挿しない補間器で定義域の外を問い合わせたとき、および NaN を問い合わせたときは `None` です。埋める値は呼び出し時に渡します: `eval_or(point, fill)`、または `eval_batch_into(points, out, fill)`。
- **定義域の外での振る舞いは、外挿のモード `Extrapolate` で決まります。** `with_extrapolate` で設定します。既定は SciPy に合わせてあります。
  - `Extend`（端の区間の式をそのまま延ばす）: `CubicSpline`、`PchipInterpolator`、`CubicHermiteSpline`、`BSpline1D`、`BSplineCurve`、`PPoly`、`BSpline`。周期境界条件で作ったスプラインは `Periodic` です。
  - `Forbid`（値なし）: `Akima1D`、`Nearest1D`、`RegularGridInterpolator`。
  - `Clamp`（端の値）: `Linear1D`。`numpy.interp` と同じです。
- **入力は検査され、SciPy が別のものを計算する入力は拒否されます。** 1 変数では座標が狭義単調増加である必要があります（SciPy の `interp1d` は並べ替えます）。データはどこでも有限である必要があり、NaN や無限大を含むデータはエラー（`InterpolateError::NonFinite`）になります。データの NaN が補間の値に紛れ込むことはありません。
- **規則格子の `cubic` と `quintic` は厳密に解きます。** SciPy の既定は反復解法で解くため、値がスプラインからずれます。ずれは、滑らかなデータで最大値の 1e-5〜1e-2（2 次元・3 次元の格子、等間隔でも不等間隔でも）、乱数の値では値そのものの大きさに達します。interpx はスプラインそのものを返し、SciPy の `method="cubic_legacy"` と丸め誤差の範囲で一致します。
- **`RBFInterpolator` に既定のカーネルはありません。** 系は LU 分解で解きます（SciPy は対称な分解。値の差は 1e-12 未満）。重複した点と、カーネルが求める最小の次数より低い多項式の次数は拒否します（次数については、SciPy は警告を出すだけです）。
- **`LinearNDInterpolator` は、同じ円の上にある 4 点以上（球の上なら 5 点以上、…）が作る区画を、独自の規則で分割します。** SciPy の三角形分割（Qhull）は別の分け方をするので、そのような区画の中では値が SciPy と違います。格子の上の点がこれに当たります。座標が同じ点が複数あるときは、最初に渡された点の値を使います（Qhull はそのうちの 1 つを残しますが、最初の点とは限りません）。また、SciPy が受け付けない 1 次元の点も扱えます。
- **値は実数です。** 複素数は扱いません。値を固定長の配列 `[T; M]` にすると、`M` 列をまとめて補間できます。
- **補間の方式は構築時に決めます。** 評価のたびには変えられません。

## 例

```rust
use interpx::{BoundaryCondition, CubicSpline, Extrapolate, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.5, 5.0];
    let y = [1.0, 2.7, 5.8, 4.0, 2.2];

    // SciPy の `CubicSpline(x, y)`。
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

- **補間器はすべて自分のデータを所有し**、`Send + Sync + 'static` です。保持しておくことも、別のスレッドへ渡すことも、共有することもできます。コンストラクタは渡されたものを複製します。データの一部を渡されたまま保持する型（`PPoly`、`BSpline`、`RegularGridInterpolator`、`NearestNDInterpolator`、`LinearNDInterpolator`、`RBFInterpolator`）には `from_owned` もあり、そのデータを `Vec` で受け取って複製せずに保持します。
- **座標は `f32` か `f64` です**（1 変数、曲線、規則格子）。散布点は `f64` だけです。
- **規則格子の値は C 順です。** 最後の軸が最も速く変わります。NumPy が `values[i, j, k]` を並べる順と同じです。
- **NaN と無限大**: データに含まれるとき、構築時のエラーになります。NaN を問い合わせると値はありません（`None`、または埋める値）。無限大の問い合わせは定義域の外として扱います。`Forbid` と `Periodic` では値なし、`Clamp` では端の値です。`Extend` では端の区間の式をその点で評価するので、値が無限大や NaN になることがあります。
- **全点の `RBFInterpolator` は、構築の間 `N` × `N` の行列を持ちます。** 10 838 点で 940 MB です。近傍の点だけを使う設定（`RbfOptions::with_neighbors`）では、代わりに問い合わせごとに小さな系を解きます。`reciprocal_condition` は、全点の系がどれだけ正確に解けたかを示します。
- **評価はメモリを確保しません。** 例外は、`Vec` を返す `eval_batch` と `eval_batch_par`、および近傍を使う `RBFInterpolator` のスレッドごとの最初の評価（そのスレッドが持ち続ける作業領域を確保します）です。
- **0.x のリリース**: マイナーバージョン（0.1 から 0.2）で公開 API が変わることがあります。ここでいう公開 API は、既定のフィーチャーで公開されるものを指します。
- **この crate にないもの**: `interp2d`（SciPy から削除済み）と、型としての `interp1d`（その種類は `Linear1D`、`Nearest1D`、`BSpline1D` で扱えます）。関数の `griddata` もありません。代わりに、`method` に当たる補間器を作って評価してください（`nearest` は `NearestNDInterpolator`、`linear` は `LinearNDInterpolator`、`cubic` は `CloughTocher2DInterpolator`。最後の 1 つはまだありません）。まだないもの: `CloughTocher2DInterpolator`、規則格子と散布点の補間の導関数、`roots`、`NdPPoly`、`BarycentricInterpolator` と `KroghInterpolator`、FITPACK 系（`UnivariateSpline`、`splrep`、`RectBivariateSpline`）、実行時に決める次元と出力の幅。

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

### エラー: `InterpolateError`、`AxisFault`、`ParamReason`、`SolverError`

コンストラクタは `InterpolateError` を返し、どの入力のどこが悪いかを示します。`InsufficientPoints`、`LengthMismatch`、`NonFinite`、`NotSorted`、`DuplicatePoints`、`NotPeriodic`、`InvalidParameter`（`ParamReason` つき）、`SolverFailed`（`SolverError` つき）、`AllocationFailed`、`GridAxis`（規則格子の 1 本の軸について。`AxisFault` つき）があります。

```rust
use interpx::{
    AxisFault, BSpline1D, BoundaryCondition, GridMethod, InterpolateError, Linear1D,
    ParamReason, RBFInterpolator, RbfKernel, RegularGridInterpolator, SolverError,
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

区分多項式（SciPy の `PPoly`）です。区切りの点と、区間ごとの係数（次数の高い順）から作ります。`BSpline` と同じく、`derivative`、`antiderivative`、`eval_derivative`、`integrate` を持ちます。

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

ノットと制御点から作る B スプライン（SciPy の `BSpline`）です。次数は `BSPLINE_MAX_DEGREE`（15）までです。`PPoly` と相互に変換できます。

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

線形補間（`numpy.interp`、`interp1d(kind="linear")`）です。データ点での値は、データの値と厳密に一致します。

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

近くのデータ点の値を返します（`interp1d` の `kind` が `"nearest"`、`"nearest-up"`、`"previous"`、`"next"` のもの）。

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

C² 級の 3 次スプライン（SciPy の `CubicSpline`）です。`BoundaryCondition` は `NotAKnot`、`Periodic`、または両端に `EndCondition`（`NotAKnot`、`FirstDerivative`、`SecondDerivative`）を指定する `Ends` です。`natural()` と `clamped()` は、SciPy の `"natural"` と `"clamped"` に当たります。

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

形を保つ C¹ 級の 3 次補間（SciPy の `PchipInterpolator`）です。データが単調な区間では補間も単調で、データの値の範囲を越えません。

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

秋間の C¹ 級の 3 次補間（SciPy の `Akima1DInterpolator`、`method="akima"` または `"makima"`）です。各点の傾きを、近くの 5 点から決めます。

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

与えた値と傾きを通る 3 次補間（SciPy の `CubicHermiteSpline`）です。

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

データ点を通る 1〜15 次の B スプライン（SciPy の `make_interp_spline`）です。`as_bspline` で、中身の `BSpline` を参照できます。

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

`D` 次元の点列を通るパラメトリック曲線です。`new` は、決め方（`Uniform`、`ChordLength`、`Centripetal`）に従ってパラメータを求めます。`with_parameters` は、`make_interp_spline(t, points)` と同じくパラメータを受け取ります。

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

`D` 次元の直交格子の上の値を補間します（SciPy の `RegularGridInterpolator`）。方式は `Nearest`、`Linear`、`Cubic`、`Quintic`、`Pchip` です。軸の間隔は不均一でもよく、降順で与えることもできます。

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

`D` 次元で、いちばん近いデータ点の値を返します（SciPy の `NearestNDInterpolator`）。`rescaled` は、各座標を同じ範囲にそろえてから距離を測ります（`rescale=True`）。`with_distance_bound` は、ある距離より遠い問い合わせを値なしにします。

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

`D` 次元で、点の Delaunay 分割の上で線形に補間します（SciPy の `LinearNDInterpolator`）。点の凸包の外では値がありません。`rescaled` は、各座標を同じ範囲にそろえてから三角形分割します（`rescale=True`）。`with_simplex_tolerance` は、単体の内側とみなす範囲を広げます（SciPy 1.18 の `simplex_tolerance`。SciPy では評価のたびに渡します）。

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

### `RBFInterpolator`、`RbfKernel`、`RbfDegree`、`RbfOptions`

`D` 次元の放射基底関数による補間（SciPy の `RBFInterpolator`）です。カーネルは SciPy と同じです。形状パラメータのある 4 つのカーネルでは、`epsilon` を指定します。`RbfOptions` で、平滑化、多項式の項の次数（`RbfDegree`）、近傍の点の数を設定します。

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

## 性能

SciPy 1.18 と、ほかの純粋な Rust の補間 crate との比較（時間とメモリ。マシン、バージョン、データ、コマンドつき）: [docs/performance.md](https://github.com/YUKIKEDA/interpx/blob/main/docs/performance.md)。

## ライセンス

[Apache License, Version 2.0](LICENSE-APACHE) と [MIT license](LICENSE-MIT) のどちらかを選んで利用できます。
