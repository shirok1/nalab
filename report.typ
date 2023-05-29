#import "template.typ": *
#set text(font: "Source Han Serif SC", lang: "zh")

// Take a look at the file `template.typ` in the file panel
// to customize this template and discover how it works.
#show: project.with(
  title: "数值计算原理实习作业",
  authors: (
    (name: "秦泽钊", affiliation: "理学院 数学三班", phone: "210810311"),
  ),
)

= 多项式插值

下面代码片段中的 `self.values` 均为提前计算的节点处的函数值。

== 拉格朗日插值

使用形如 $ l_k(x) = ((x-x_0) dots.c (x-x_(k-1))(x-x_(k+1)) dots.c (x-x_n))/((x_k-x_0) dots.c (x_k-x_(k-1))(x_k-x_(k+1)) dots.c (x_k-x_n)) $ 的插值基函数的插值多项式即为拉格朗日插值，即 $ L_k(x) = sum_(k=0)^n y_k l_k(x) = sum_(k=0)^n y_k omega(x)/((x-x_k) omega'_k), \ omega(x)=(x-x_0)(x-x_1) dots.c (x-x_n), \ omega'_k = (x_k-x_0)(x_k-x_1) dots.c (x_k-x_(k-1))(x_k-x_(k+1)) dots.c (x_k-x_n) $

注意到 $omega'_k$ 是一个与 $x$ 无关，仅与插值节点有关的值，因此可以在插值初始化时计算完毕。

```rs
let values = nodes.iter().copied().map(func).collect();

let omega_prime = nodes.iter().enumerate().map(|(i, x_i)| {
    nodes.iter().take(i).chain(nodes.iter().skip(i + 1))
        .fold(1.0, |acc, x_j| { acc * (x_i - x_j) })
}).collect::<Vec<f64>>();
```

另外，在计算 $L_k(x)$ 时，如果 $x$ 存在于节点中（即存在 $x=x_k$），$omega(x)$ 和 $x-x_k$ 会同时为 $0$，此时进行除法就会产生错误，因此需要进行判断并用直接计算代替。

```rs
let omega = self.nodes.iter()
    .map(|x_i| x - x_i)
    .fold(1.0, |acc, x_i| acc * x_i);

let result = self.nodes.iter().enumerate()
    .map(|(k, x_k)| {
        let y_k = self.values[k];
        if omega.is_normal() && self.omega_prime[k].is_normal() {
            y_k * omega / ((x - x_k) * self.omega_prime[k])
        } else {
            // Can't use cached omega_prime
            let x_k = self.nodes[k];
            let l_k = self.nodes[..k].iter().chain(&self.nodes[k + 1..])
                .fold(1.0, |acc, x_i|
                    acc * (x - x_i) / (x_k - x_i));
            y_k * l_k
        }
    }).sum::<f64>();
```

== 切比雪夫多项式零点

直接编写 `chebyshev_zeros` 计算 $-cos((2k + 1)/(2n)pi)$ 即可，插值代码中不应当对插值节点作顺序以外的假设。

== 分段线性插值

找到离 $x$ 最近的一或两个节点，进行线性插值即可。其中找最近点可以使用二分查找等加速方法（代码中未实现）。

```rs
let (a, b) = find_between(self.nodes, x);

let x_a = self.nodes[a];
let x_b = self.nodes[b];
let y_a = self.values[a];
let y_b = self.values[b];

let result = y_a + (y_b - y_a) * (x - x_a) / (x_b - x_a);
```

== 三次样条插值

由题意可知，题目符合 “已知两端的一阶导数值” 的定义，因此要求插值多项式，就要先求出以下非齐次线性方程组的解：$ mat(
  2, lambda_0;
  mu_1, 2, lambda_1;
  "", dots.down, dots.down, dots.down;
  "", "", mu_(n-1), 2, lambda_(n-1);
  "", "", "", mu_n, 2;
) mat(M_0; M_1; dots.v; M_(n-1); M_n) =  mat(d_0; d_1; dots.v; d_(n-1); d_n) $

计算得到 $b$ 之后即可对每一个传入的 $x$ 找到应该出现在的分段，并完成插值的计算：$ x = 1/(6h) dot (b_a (x_b - x)^3 + b_b (x - x_a)^3) + (y_a / h - (b_a h) / 6) (x_b - x) + (y_b / h - (b_b h) / 6) (x - x_a) $

遗憾的是，因为大小无法在编译期确定，代码中的解线性方程组无法复用迭代求解线性方程组中的代码。

== 图像展示

#figure(
  image("interpolate.svg"),
  caption: [
    原函数、插值函数、拟合函数的图像展示（正半轴部分）
  ],
)


== 效果检验

#figure(
  table(
    columns: (auto, 1fr),
    align: center,
    [插值方法], [最大误差（与实际值）],
    "拉格朗日插值", $9.0081 times 10^(0)$,
    "拉格朗日插值，切比雪夫多项式原点", $5.9993 times 10^(-3)$,
    "线性插值", $5.9993 times 10^(-3)$,
    "三次样条插值", $2.2762 times 10^(-5)$,
),
  caption: [
    多项式插值结果表
  ],
)

= 数值积分

== 复合梯形公式、复合辛普森公式

按照公式实现即可。代码中分别实现为 `nalab::task2::composite_trapezoid` 和 `nalab::task2::composite_simpson`

== 龙贝格积分

在实际实现时，利用了 Rust 语言的 “切片引用” 特性 #footnote[https://doc.rust-lang.org/std/primitive.slice.html] 设计了包含依赖类型思想的缓存结构。众所周知，在计算龙贝格积分中的值 $T_m^((k)) (m!=0)$ 时会利用 $T_(m-1)^((k))$ 和 $T_(m-1)^((k+1))$ 两个值。按照递推关系可以画出如下的 $T$ 表：

#figure(
  table(
    columns: (auto, auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
    align: center,
    [$k$], [$h$], [$T_0^((k))$], [$T_1^((k))$], [$T_2^((k))$], [$T_3^((k))$], [$T_4^((k))$], [$dots.c$],
    $0$, $b-a$, $T_0^((0))$, [], [], [], [] , $dots.c$,
    $1$, $(b-a)/2$, $T_0^((1))$, $T_1^((0))$, [], [], [], $dots.c$,
    $2$, $(b-a)/4$, $T_0^((2))$, $T_1^((1))$, $T_2^((0))$, [], [], $dots.c$,
    $3$, $(b-a)/8$, $T_0^((3))$, $T_1^((2))$, $T_2^((1))$, $T_3^((0))$, [], $dots.c$,
    $4$, $(b-a)/16$, $T_0^((4))$, $T_1^((3))$, $T_2^((2))$, $T_3^((1))$, $T_4^((0))$, $dots.c$,
    $dots.v$, $dots.v$, $dots.v$, $dots.v$, $dots.v$, $dots.v$, $dots.v$, $dots.down$,
),
  caption: [
    $T$ 表
  ],
)

从表中可以看出，如果直接使用递归的写法，会导致有大量的值被重复计算。假如使用一般的函数记忆化方法，会产生函数难以进行散列的问题。此时我们可以看出一个现象：计算一个值仅需要其前置项（实际上，这个表按照递推的计算顺序来排序）而不需要任何后置项。因此，在递归计算时，我们可以不断“缩小”所需的切片的范围。对于 $T_(m-1)^((k+1))$，所需的切片刚好是除去最后一项（即当前项）的剩余部分，而 $T_(m-1)^((k))$ 需要的切片更小，可以往上排除一整行：

```rs
fn romberg_rec(f: fn(f64) -> f64, zone: Range<f64>, partition: usize, partition_exp: usize, acceleration: usize, buffer: &mut [Option<f64>]) -> f64 {
    let line = partition_exp + acceleration;
    assert_eq!(buffer.len(), (line + 1) * (line + 2) / 2 - partition_exp);
    let (current, rest) = buffer.split_last_mut().unwrap();
    *(current.get_or_insert_with(|| {
        if acceleration == 0 {
            composite_trapezoid(f, zone, partition * 2usize.pow(partition_exp as u32))
        } else {
            let numerator = 4u32.pow(acceleration as u32);
            let denominator = numerator - 1;
            let rest_len = rest.len();
            ((numerator as f64) * romberg_rec(f, zone.clone(), partition, partition_exp + 1, acceleration - 1, rest)
                - romberg_rec(f, zone.clone(), partition, partition_exp, acceleration - 1, rest[0..rest_len - line].as_mut()))
                / (denominator as f64)
        }
    }))
}
```

不过，虽然上面的写法非常符合函数式编程的要求，时间复杂度也不高，空间复杂度却达到了 $O(k^2)$，原因是递归调用的模式下很难及时的“丢弃”不需要的空间。因此，这里另外写了一版使用递推的版本，不得不使用了 _foreach_ 循环，仅保留当前行和上一行，将空间复杂度降至 $O(k)$：

```rs

fn romberg_2(
    f: fn(f64) -> f64,
    zone: Range<f64>,
    partition: usize,
    partition_exp: usize,
    last_line: &[f64],
    mut this_line: &mut [f64],
) {
    assert_eq!(last_line.len() + 1, this_line.len());
    this_line[0] = composite_trapezoid(f, zone, partition * 2usize.pow(partition_exp as u32));
    for (i, (t_m_minus_1_k, [t_m_minus_1_k_plus_1, t_m_k])) in
    last_line.iter().zip(Cell::from_mut(this_line).as_slice_of_cells().array_windows()).enumerate() {
        let acceleration = i + 1;
        let numerator = 4u32.pow(acceleration as u32);
        let denominator = numerator - 1;
        t_m_k.set(
            ((numerator as f64) * t_m_minus_1_k_plus_1.get() - t_m_minus_1_k) / (denominator as f64)
        );
    }
}

pub fn simple_romberg_2(f: fn(f64) -> f64, zone: Range<f64>, init_partition: usize, acceleration: usize) -> f64 {
    let mut buffer_a = vec![0.; acceleration + 1];
    let mut buffer_b = vec![0.; acceleration + 1];

    let mut last_line = &mut buffer_a[..];
    let mut this_line = &mut buffer_b[..];

    for i in 0..=acceleration {
        romberg_2(f, zone.clone(), init_partition, i, &last_line[0..i], &mut this_line[0..i + 1]);
        std::mem::swap(&mut last_line, &mut this_line);
    }
    last_line[acceleration]
}
```

== 复合高斯公式

本程序中将其实现为先将区间分为若干个子区间，将他们分别缩放到 $(-1, 1)$，然后再使用高斯——勒让德求积公式中可以预先算出的节点与权重完成数值积分，即实际上实现了 $ integral_a^b f(x) rho(x) diff x = integral_a^b q(x) rho(x) diff x = sum_(k=0)^n A_k f(x_k) $



为了提高精度，这里的节点与权重没有使用手动输入小数值，而是使用了编译期常值计算。

== 效果检验

#figure(
  table(
    columns: (auto, 1fr),
    align: center,
    [数值积分方法], [误差（与实际值）],
    "复合梯形公式，等分 40 份", $9.0081 times 10^(0)$,
    "复合辛普森公式，等分 20 份", $5.9993 times 10^(-3)$,
    "复合辛普森公式，等分 20 份，龙贝格加速一次", $5.9993 times 10^(-3)$,
    "复合辛普森公式，等分 10 份，龙贝格加速两次", $2.2762 times 10^(-5)$,
    "复合辛普森公式，等分 5 份，龙贝格加速三次", $3.5838 times 10^(-7)$,
    "复合的 2 点高斯公式，等分 20 份", $3.9992 times 10^(-3)$,
    "复合的 4 点高斯公式，等分 10 份", $1.8190 times 10^(-12)$,
),
  caption: [
    数值积分结果表
  ],
)

可以看出，这里效果最好的数值积分方法是 “复合的 4 点高斯公式，等分 10 份” 的场景。值得一提的是，无论选择表里的哪种数值积分方法，都只会实际调用 $f$ 四十次。

= 迭代求解线性方程组

此处利用了 Rust 的泛型系统 #footnote[https://doc.rust-lang.org/book/ch10-01-syntax.html] 来自动在编译期适配不同大小的矩阵。另外，本次需要实现的迭代法都没有中间状态，因此可以通过 Rust 的特性系统抽象出不同的求解方法通用的迭代执行器 `nalab::task3::abstract_iteration`，使得不同的求解方法在代码上几乎只有读入旧值、返回新值的闭包的区别。

== 雅可比迭代法、高斯——赛德尔迭代法

这两种方法的的实现都较为简单，直接按照公式：$ x_("Jacobi")^((k+1)) = D^(-1) ((L + U) x^((k)) + b) \ x_("GS")^((k+1)) = (D - L)^(-1) (U x^((k)) + b) $ 进行迭代即可。

== 超松弛（SOR）迭代法

需要先手动求出该矩阵的谱半径 $rho(A)$，然后通过 `nalab::task3::omega_opt` 按公式求出对应的最优松弛因子 $omega_("opt")$，接下来的流程与高斯——赛德尔迭代法相似，只是在迭代公式中加入了松弛因子：$ x_("SOR")^((k+1)) = (D - omega_("opt") L)^(-1) ((1 - omega_("opt")) D + omega_("opt") U) x^((k)) + omega_("opt") (D - omega_("opt") L)^(-1) b $

== 最速梯度下降法

首先计算出剩余向量 $r^(k) = b - A x^((k))$，然后按照公式 $ x_("SD")^((k+1)) = x^((k)) + alpha^(k) r^(k) $ 进行迭代，其中 $alpha^(k) = (r^(k), r^(k)) / (A r^(k), r^(k))$。

== 效果检验

#figure(
  table(
    columns: (auto, 1fr, 1fr, 1fr),
    align: center,
    [方法], [步数], [残差（剩余向量）], [误差（与实际值）],
    "雅可比迭代法", $45$, $9.9246 times 10^(-9)$, $3.9698 times 10^(-9)$,
    "高斯——赛德尔迭代法", $13$, $3.2293 times 10^(-9)$, $6.9968 times 10^(-10)$,
    "超松弛迭代法", $12$, $5.9476 times 10^(-9)$, $1.0519 times 10^(-9)$,
    "最速梯度下降法", $22$, $9.1539 times 10^(-9)$, $3.7915 times 10^(-9)$,
),
  caption: [
    迭代求解线性方程组结果表
  ],
)

因为结束条件本身设置成了残差，因此最终的残差都相差无几，但此时 “高斯——赛德尔迭代法” 同时拥有极低的步数和最低的误差，可以认为是最好的选择。
