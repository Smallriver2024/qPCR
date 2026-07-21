// Type declarations for the `jstat` package (no official @types/jstat exists).
// Covers the statistical APIs typically needed by this project; the index
// signature keeps the rest of the library accessible without full typings.
declare module "jstat" {
  interface ContinuousDistribution {
    pdf(x: number, ...params: number[]): number;
    cdf(x: number, ...params: number[]): number;
    inv(p: number, ...params: number[]): number;
    mean(...params: number[]): number;
    median(...params: number[]): number;
    mode(...params: number[]): number;
    variance(...params: number[]): number;
    sample(...params: number[]): number;
  }

  interface JStat {
    // Descriptive statistics
    mean(arr: number[] | number[][]): number;
    median(arr: number[] | number[][]): number;
    mode(arr: number[] | number[][]): number;
    sum(arr: number[] | number[][]): number;
    sumsqrd(arr: number[]): number;
    sumsqerr(arr: number[]): number;
    sumrow(arr: number[][]): number;
    product(arr: number[]): number;
    min(arr: number[] | number[][]): number;
    max(arr: number[] | number[][]): number;
    range(arr: number[]): number;
    variance(arr: number[] | number[][], flag?: boolean): number;
    pooledvariance(arr: number[][]): number;
    deviation(arr: number[]): number[];
    stdev(arr: number[] | number[][], flag?: boolean): number;
    pooledstdev(arr: number[][]): number;
    meandev(arr: number[]): number;
    meddev(arr: number[]): number;
    skewness(arr: number[]): number;
    kurtosis(arr: number[]): number;
    coeffvar(arr: number[]): number;
    quartiles(arr: number[]): number[];
    quantiles(arr: number[], quantiles: number[]): number[];
    percentile(arr: number[], k: number): number;
    percentileOfScore(arr: number[], score: number, kind?: string): number;
    covariance(arr1: number[], arr2: number[]): number;
    corrcoeff(arr1: number[], arr2: number[]): number;
    spearmancoeff(arr1: number[], arr2: number[]): number;
    geomean(arr: number[]): number;
    cumsum(arr: number[]): number[];
    cumprod(arr: number[]): number[];
    diff(arr: number[]): number[];
    rank(arr: number[]): number[];

    // Hypothesis tests
    /** One-sample / two-sample t-test score */
    tscore(...args: number[]): number;
    /** Two-sided p-value for a t score with given degrees of freedom */
    ttest(...args: number[]): number;
    /** One-way ANOVA F score */
    anovafscore(...args: (number[] | number[][])[]): number;
    /** p-value for an ANOVA F score */
    anovaftest(...args: number[][][]): number;
    /** Two-sample F test for equal variances */
    ftest(fscore: number, df1: number, df2: number): number;
    /** Tukey's range test (single pairwise comparison) */
    tukeyhsd(arrays: number[][]): number;
    qscore(...args: number[]): number;
    qtest(...args: number[]): number;

    // Special functions
    gammap(a: number, x: number): number;
    gammaln(x: number): number;
    gammafn(x: number): number;
    betafn(a: number, b: number): number;
    betaln(a: number, b: number): number;
    ibetainv(p: number, a: number, b: number): number;
    ibeta(x: number, a: number, b: number): number;
    factorial(n: number): number;
    combination(n: number, k: number): number;
    permutation(n: number, k: number): number;
    randn(n?: number, m?: number): number | number[] | number[][];
    randg(shape: number, n?: number, m?: number): number | number[] | number[][];

    // Distributions
    normal: ContinuousDistribution;
    studentt: ContinuousDistribution;
    centralF: ContinuousDistribution;
    chisquare: ContinuousDistribution;
    exponential: ContinuousDistribution;
    beta: ContinuousDistribution;
    gamma: ContinuousDistribution;
    weibull: ContinuousDistribution;
    cauchy: ContinuousDistribution;
    lognormal: ContinuousDistribution;
    pareto: ContinuousDistribution;
    uniform: ContinuousDistribution;
    triangular: ContinuousDistribution;
    invgamma: ContinuousDistribution;
    kumaraswamy: ContinuousDistribution;
    laplace: ContinuousDistribution;
    noncentralt: ContinuousDistribution;
    binomial: ContinuousDistribution & { pdf(k: number, n: number, p: number): number };
    negbin: ContinuousDistribution;
    hypgeom: ContinuousDistribution;
    poisson: ContinuousDistribution;

    [key: string]: unknown;
  }

  const jStat: JStat;
  export = jStat;
}
