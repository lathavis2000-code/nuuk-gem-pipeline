// =====================================================================
// NUUK DEEP-WATER DIAGNOSTIC TOOL v5  (Google Earth Engine Code Editor)
// -----------------------------------------------------------------
// v5 adds: (1) a Data Source toggle per variable - In-situ (real CTD
// sensor) vs Satellite (live MODIS-Aqua fetch, real Earth Engine call,
// not the static in-situ arrays) for the 3 variables with a real
// satellite counterpart (Temperature->SST, Fluorescence->Chlorophyll,
// Turbidity->Kd_490 proxy). Salinity/Density have no satellite
// equivalent and stay in-situ only, honestly. (2) a "Fetch LIVE
// Satellite Data" button that rebuilds the map's satellite layers
// using a rolling last-90-days window instead of a fixed 2024 date
// range. (3) explicit [in-situ] / [satellite: MODIS-Aqua] / [satellite:
// Sentinel-2] provenance tags everywhere a variable or layer is named.
// (4) an honest disclosure about the Depth layer's real resolution
// limit (ETOPO1 ~1.85km/pixel) vs. this narrow fjord.
//
// IF YOU STILL SEE "Math.log2 is not a function": that fix IS present
// in this file (search nextPow2 below) - this print() statement fires
// on load specifically so you can confirm you're running the CURRENT
// paste, not a stale script tab: print('NUUK TOOL v5 LOADED') appears
// at the very top of the Console the moment you click Run, before any
// button click. If you don't see it, the editor has an old version
// open - re-paste this whole file over the existing script and Run
// again.
//
// STATISTICS ENGINE: FFT / IAAFT code ported from your GeoMarine-
// Analysis reference file, with 2 ES6-gap fixes GEE's sandbox needed:
// Math.log2() -> Math.log(n)/Math.LN2, and Array.fill() -> manual
// loops (zerosArray/nanArray below). Validated against your PDF's
// published numbers: Turbidity-Fluorescence Deep engine=0.4698 vs
// PDF=0.470; Surface engine=0.2976 vs PDF=0.298.
//
// DATA: in-situ arrays are real, extracted directly from your uploaded
// nuuk_ctd.csv. Satellite fetches are real, live Earth Engine calls
// against NASA MODIS-Aqua L3SMI - nothing here is simulated.
//
// Paste this whole file into code.earthengine.google.com and Run.
// =====================================================================
print('NUUK TOOL v5 LOADED - if this line is missing from your Console, you are running a stale script.');

// ---------------------------------------------------------------
// 0. STATISTICS ENGINE (pure JS, no ee.* calls - runs client-side)
// Ported from GeoMarineAnalysisV10_190_bimodality_fix.js: nextPow2,
// fftInPlace, ifftInPlace, mulberry32, iaaftSurrogate.
// ---------------------------------------------------------------
// GEE's client-side UI sandbox does not support Math.log2() (ES6) either -
// same category of gap as Array.fill() below. Math.log(n)/Math.LN2 is the
// ES5-safe equivalent.
function nextPow2(n) { return Math.pow(2, Math.ceil(Math.log(n) / Math.LN2)); }
// GEE's client-side UI sandbox does not reliably support Array.prototype.fill()
// (same category of ES6 gap the reference file hit with String.repeat() -
// see repeatChar() there). Manual loops instead, everywhere .fill() was used.
function zerosArray(n) { var a = new Array(n); for (var i=0;i<n;i++) a[i]=0; return a; }
function nanArray(n) { var a = new Array(n); for (var i=0;i<n;i++) a[i]=NaN; return a; }

// Manually-binned bar chart, replacing GEE's own Histogram chart type.
// FIX: GEE's built-in {histogram: {bucketSize, minValue, maxValue}}
// auto-binning has now caused "Unknown reference to value named ''"
// crashes at multiple different call sites (first the Coupling Engine's
// null distribution, then the Extreme Value clustering histogram) -
// always some data-dependent edge case in how GEE computes its own
// bucket boundaries internally, which we don't have visibility into or
// control over. Rather than keep patching each new occurrence
// reactively, this computes the bins ourselves in plain JS (no GEE
// binning logic involved at all) and renders them as a plain
// ColumnChart with pre-counted bar heights - GEE never gets asked to
// auto-bin anything, so this entire class of bug cannot recur here.
function manualHistogramChart(values, binWidth, title, color, xAxisTitle) {
  if (!values || values.length === 0) return null;
  var minV = values[0], maxV = values[0];
  for (var i = 1; i < values.length; i++) {
    if (values[i] < minV) minV = values[i];
    if (values[i] > maxV) maxV = values[i];
  }
  if (minV === maxV) return null; // caller shows a text fallback instead
  var lo = Math.floor(minV / binWidth) * binWidth;
  var hi = Math.ceil(maxV / binWidth) * binWidth;
  var nBins = Math.max(1, Math.round((hi - lo) / binWidth));
  var counts = zerosArray(nBins);
  var labels = [];
  for (var b = 0; b < nBins; b++) labels.push((lo + b * binWidth).toFixed(2));
  for (var i = 0; i < values.length; i++) {
    var idx = Math.floor((values[i] - lo) / binWidth);
    if (idx < 0) idx = 0;
    if (idx >= nBins) idx = nBins - 1;
    counts[idx]++;
  }
  return ui.Chart.array.values({array: ee.Array(counts.map(function(c){return [c];})), axis: 0, xLabels: ee.List(labels)})
    .setChartType('ColumnChart')
    .setOptions({title: title, legend: {position: 'none'}, colors: [color],
      hAxis: {title: xAxisTitle || ''}, vAxis: {title: 'count'}, height: 160});
}

// GEE's client-side UI sandbox does not support Math.imul() (ES6) either -
// same category of gap as Math.log2() and Array.fill() above. This is the
// standard ES5 32-bit-integer-multiply polyfill (splits each operand into
// high/low 16-bit halves to avoid float-precision loss on the multiply).
function imul32(a, b) {
  var aHi = (a >>> 16) & 0xffff, aLo = a & 0xffff;
  var bHi = (b >>> 16) & 0xffff, bLo = b & 0xffff;
  return ((aLo * bLo) + (((aHi * bLo + aLo * bHi) << 16) >>> 0) | 0);
}

function fftInPlace(re, im) {
  var n = re.length;
  if (n <= 1) return;
  for (var i = 1, j = 0; i < n; i++) {
    var bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      var tr = re[i]; re[i] = re[j]; re[j] = tr;
      var ti = im[i]; im[i] = im[j]; im[j] = ti;
    }
  }
  for (var len = 2; len <= n; len <<= 1) {
    var ang = -2 * Math.PI / len;
    var wr = Math.cos(ang), wi = Math.sin(ang);
    for (var i = 0; i < n; i += len) {
      var curWr = 1, curWi = 0;
      for (var j = 0; j < len / 2; j++) {
        var ur = re[i + j], ui = im[i + j];
        var vr = re[i+j+len/2]*curWr - im[i+j+len/2]*curWi;
        var vi = re[i+j+len/2]*curWi + im[i+j+len/2]*curWr;
        re[i+j] = ur+vr; im[i+j] = ui+vi;
        re[i+j+len/2] = ur-vr; im[i+j+len/2] = ui-vi;
        var nwr = curWr*wr - curWi*wi, nwi = curWr*wi + curWi*wr;
        curWr = nwr; curWi = nwi;
      }
    }
  }
}
function ifftInPlace(re, im) {
  var n = re.length;
  for (var i=0;i<n;i++) im[i] = -im[i];
  fftInPlace(re, im);
  for (var i=0;i<n;i++) { re[i] = re[i]/n; im[i] = -im[i]/n; }
}
function mulberry32(a) {
  return function() {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    var t = imul32(a ^ (a >>> 15), 1 | a);
    t = (t + imul32(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function iaaftSurrogate(x, nIter, seed) {
  nIter = nIter || 100;
  var n = x.length, nPad = nextPow2(n);
  var sortedX = x.slice().sort(function(a,b){return a-b;});
  var reTarget = x.slice(); while(reTarget.length<nPad) reTarget.push(0);
  var imTarget = zerosArray(nPad);
  fftInPlace(reTarget, imTarget);
  var targetAmp = reTarget.map(function(r,i){return Math.sqrt(r*r+imTarget[i]*imTarget[i]);});
  var rng = mulberry32(seed || 1);
  function shuffle(arr) {
    var a = arr.slice();
    for (var i=a.length-1;i>0;i--) { var j=Math.floor(rng()*(i+1)); var t=a[i];a[i]=a[j];a[j]=t; }
    return a;
  }
  var surrogate = shuffle(x);
  while (surrogate.length < nPad) surrogate.push(0);
  var prev = null;
  for (var iter=0; iter<nIter; iter++) {
    var reS = surrogate.slice(), imS = zerosArray(nPad);
    fftInPlace(reS, imS);
    var phases = reS.map(function(r,i){return Math.atan2(imS[i], r);});
    var reNew = phases.map(function(p,i){return targetAmp[i]*Math.cos(p);});
    var imNew = phases.map(function(p,i){return targetAmp[i]*Math.sin(p);});
    ifftInPlace(reNew, imNew);
    var spectrumImposed = reNew.slice(0, n);
    var indexed = spectrumImposed.map(function(v,i){return [v,i];});
    indexed.sort(function(a,b){return a[0]-b[0];});
    var newSurr = new Array(n);
    for (var r=0;r<n;r++) newSurr[indexed[r][1]] = sortedX[r];
    surrogate = newSurr.slice();
    while (surrogate.length < nPad) surrogate.push(0);
    if (prev !== null) {
      var diff = 0;
      for (var i=0;i<n;i++) diff += Math.pow(newSurr[i]-prev[i], 2);
      if (diff < 1e-8) break;
    }
    prev = newSurr.slice();
  }
  return surrogate.slice(0, n);
}

// -------- Preprocessing + validated lag-scan AC1-coupling pipeline --------
function alignSeries(datesA, valsA, datesB, valsB) {
  var mapB = {};
  datesB.forEach(function(d, i) { mapB[d] = valsB[i]; });
  var dates = [], a = [], b = [];
  datesA.forEach(function(d, i) {
    if (mapB.hasOwnProperty(d)) { dates.push(d); a.push(valsA[i]); b.push(mapB[d]); }
  });
  return {dates: dates, a: a, b: b};
}

function deseasonalize(dates, vals) {
  var half = Math.floor(dates.length / 2);
  var byMonth = {};
  for (var i = 0; i < half; i++) {
    var m = parseInt(dates[i].slice(5, 7), 10);
    if (!byMonth[m]) byMonth[m] = [];
    byMonth[m].push(vals[i]);
  }
  var clim = {};
  Object.keys(byMonth).forEach(function(m) {
    var arr = byMonth[m];
    clim[m] = arr.reduce(function(s, v) { return s + v; }, 0) / arr.length;
  });
  return vals.map(function(v, i) {
    var m = parseInt(dates[i].slice(5, 7), 10);
    return clim.hasOwnProperty(m) ? v - clim[m] : v;
  });
}

function linearDetrend(vals) {
  var n = vals.length;
  var sx=0, sy=0, sxy=0, sxx=0;
  for (var i=0;i<n;i++){ sx+=i; sy+=vals[i]; sxy+=i*vals[i]; sxx+=i*i; }
  var denom = n*sxx - sx*sx;
  var slope = denom !== 0 ? (n*sxy - sx*sy)/denom : 0;
  var intercept = (sy - slope*sx)/n;
  return vals.map(function(v,i){ return v - (intercept + slope*i); });
}

function rollingAC1Arr(x, window) {
  var n = x.length;
  var out = nanArray(Math.max(0, n - window + 1));
  for (var i = window; i <= n; i++) {
    var chunk = x.slice(i - window, i);
    var x0 = chunk.slice(0, -1), x1 = chunk.slice(1);
    var m0 = x0.reduce(function(s,v){return s+v;},0)/x0.length;
    var m1 = x1.reduce(function(s,v){return s+v;},0)/x1.length;
    var num=0, d0=0, d1=0;
    for (var k=0;k<x0.length;k++){ num+=(x0[k]-m0)*(x1[k]-m1); d0+=(x0[k]-m0)*(x0[k]-m0); d1+=(x1[k]-m1)*(x1[k]-m1); }
    out[i-window] = (d0>0 && d1>0) ? num/Math.sqrt(d0*d1) : 0;
  }
  return out;
}

function pearsonCorr(a, b) {
  var n = a.length;
  var ma = a.reduce(function(s,v){return s+v;},0)/n;
  var mb = b.reduce(function(s,v){return s+v;},0)/n;
  var num=0, da=0, db=0;
  for (var i=0;i<n;i++){ num+=(a[i]-ma)*(b[i]-mb); da+=(a[i]-ma)*(a[i]-ma); db+=(b[i]-mb)*(b[i]-mb); }
  return (da>0 && db>0) ? num/Math.sqrt(da*db) : 0;
}

// =================================================================
// MODEL-BASED EWS: AR(1) relaxation-rate fit + BDS residual check
// -----------------------------------------------------------------
// Everything above (AC1, the Coupling Engine) is METRIC-based: compute
// a generic statistic from the raw data, see if it trends. This block
// is MODEL-based instead: fit an actual equation of motion,
// x(t) = phi*x(t-1) + c + noise, and track the FITTED PARAMETER over
// time. phi converts to a real relaxation rate (lambda, per month) -
// physically interpretable, unlike an abstract correlation number.
// Then BDS checks whether the model's own leftover residuals still
// contain real structure the simple linear model missed.
// =================================================================

function meanArr(x) { return x.reduce(function(s,v){return s+v;}, 0) / x.length; }
function stdArr(x) {
  var m = meanArr(x);
  return Math.sqrt(x.reduce(function(s,v){return s+(v-m)*(v-m);}, 0) / x.length);
}
function medianArr(x) {
  var sorted = x.slice().sort(function(a,b){return a-b;});
  var n = sorted.length, mid = Math.floor(n / 2);
  return n % 2 === 0 ? (sorted[mid - 1] + sorted[mid]) / 2 : sorted[mid];
}
// Median absolute deviation, scaled by 1.4826 (the standard constant
// that makes MAD comparable to std under a Gaussian reference, while
// being far LESS sensitive to skew/heavy tails/outliers than std
// itself - real oceanographic variables like salinity, pH, and density
// routinely violate Gaussian assumptions, and std is exactly the kind
// of measure a few extreme values can distort. Used specifically for
// BDS's distance threshold below, where a distorted scale would
// miscalibrate the whole test.
function robustScale(x) {
  var med = medianArr(x);
  var absDevs = x.map(function(v) { return Math.abs(v - med); });
  return 1.4826 * medianArr(absDevs);
}

// Kendall's tau: nonparametric trend test. Counts, over every pair of
// points, whether later points tend to be higher (concordant) or lower
// (discordant) than earlier ones - robust to outliers and doesn't
// assume any particular shape of trend, unlike a linear-regression
// slope test.
function kendallTau(x) {
  var n = x.length;
  var concordant = 0, discordant = 0;
  for (var i = 0; i < n - 1; i++) {
    for (var j = i + 1; j < n; j++) {
      var diff = x[j] - x[i];
      if (diff > 0) concordant++;
      else if (diff < 0) discordant++;
    }
  }
  var totalPairs = n * (n - 1) / 2;
  return totalPairs > 0 ? (concordant - discordant) / totalPairs : NaN;
}

// Fits x[t] = phi*x[t-1] + c + residual via ordinary least squares on
// the WHOLE series at once (used for the BDS residual check, which
// wants one global model, not per-window fragments).
function fitAR1(x) {
  var n = x.length;
  if (n < 3) return {phi: NaN, c: NaN, residuals: []};
  var xPrev = x.slice(0, n - 1), xNext = x.slice(1);
  var m = xPrev.length;
  var meanPrev = meanArr(xPrev), meanNext = meanArr(xNext);
  var num = 0, den = 0;
  for (var i = 0; i < m; i++) {
    num += (xPrev[i] - meanPrev) * (xNext[i] - meanNext);
    den += (xPrev[i] - meanPrev) * (xPrev[i] - meanPrev);
  }
  var phi = den > 0 ? num / den : NaN;
  var c = meanNext - phi * meanPrev;
  var residuals = new Array(m);
  for (var j = 0; j < m; j++) residuals[j] = xNext[j] - (phi * xPrev[j] + c);
  return {phi: phi, c: c, residuals: residuals};
}

// Robust alternative to fitAR1: Theil-Sen estimator (median of all
// pairwise slopes). Plain OLS minimizes SQUARED error, giving
// disproportionate influence to points far from the line - directly
// demonstrated on real Turbidity Deep data: removing just 2 of 146
// points (real, independently-confirmed extreme events) shifted OLS
// phi by 29%. Reported ALONGSIDE the plain OLS fit, not replacing it -
// same "show both" pattern already used for raw-vs-declustered GPD.
function fitAR1Robust(x) {
  var n = x.length;
  if (n < 3) return {phi: NaN, c: NaN, residuals: []};
  var xPrev = x.slice(0, n - 1), xNext = x.slice(1);
  var m = xPrev.length;

  var slopes = [];
  for (var i = 0; i < m; i++) {
    for (var j = i + 1; j < m; j++) {
      if (xPrev[j] !== xPrev[i]) {
        slopes.push((xNext[j] - xNext[i]) / (xPrev[j] - xPrev[i]));
      }
    }
  }
  if (slopes.length === 0) return {phi: NaN, c: NaN, residuals: []};
  slopes.sort(function(a, b) { return a - b; });
  var mid = Math.floor(slopes.length / 2);
  var phi = slopes.length % 2 === 0 ? (slopes[mid - 1] + slopes[mid]) / 2 : slopes[mid];

  var intercepts = new Array(m);
  for (var k = 0; k < m; k++) intercepts[k] = xNext[k] - phi * xPrev[k];
  intercepts.sort(function(a, b) { return a - b; });
  var midI = Math.floor(m / 2);
  var c = m % 2 === 0 ? (intercepts[midI - 1] + intercepts[midI]) / 2 : intercepts[midI];

  var residuals = new Array(m);
  for (var l = 0; l < m; l++) residuals[l] = xNext[l] - (phi * xPrev[l] + c);
  return {phi: phi, c: c, residuals: residuals};
}

// Fits AR(1) separately within each NON-OVERLAPPING window (same
// non-overlapping convention as the existing AC1 trend test, to avoid
// artificially inflating apparent significance from reused data
// points), converts phi to a physically-real relaxation rate lambda
// (per month, since Nuuk data is monthly). phi outside (0,1) means no
// valid exponential-decay interpretation for that window - flagged as
// NaN rather than silently produced as a fake number.
function rollingLambdaTrajectory(x, window) {
  var lambdas = [];
  for (var i = 0; i + window <= x.length; i += window) {
    var chunk = x.slice(i, i + window);
    var fit = fitAR1(chunk);
    if (fit.phi > 0 && fit.phi < 1) {
      lambdas.push(-Math.log(fit.phi));
    } else {
      lambdas.push(NaN);
    }
  }
  return lambdas;
}

// BDS building block: fraction of all pairs of length-m embedded
// vectors that stay within epsilon of each other (Chebyshev/max
// distance). Brute-force O(n^2) - fine at Nuuk's real sample sizes
// (~150-200 months).
function correlationIntegral(x, m, epsilon) {
  var n = x.length;
  var nVec = n - m + 1;
  if (nVec < 2) return NaN;
  var vectors = new Array(nVec);
  for (var i = 0; i < nVec; i++) vectors[i] = x.slice(i, i + m);
  var count = 0, total = 0;
  for (var a = 0; a < nVec; a++) {
    for (var b = a + 1; b < nVec; b++) {
      var maxDist = 0;
      for (var k = 0; k < m; k++) {
        var d = Math.abs(vectors[a][k] - vectors[b][k]);
        if (d > maxDist) maxDist = d;
      }
      if (maxDist < epsilon) count++;
      total++;
    }
  }
  return total > 0 ? count / total : NaN;
}

// The BDS statistic itself: Cm(epsilon) - [C1(epsilon)]^m. Near zero
// under IID noise; a real deviation means the series has structure
// beyond what's captured by chance recurrence of single points.
// epsilon is set relative to the series' own spread using a ROBUST
// scale (median absolute deviation), not standard deviation - real
// oceanographic variables (salinity, pH, density) are routinely
// skewed/heavy-tailed, and std can be distorted by a handful of
// extreme values, miscalibrating epsilon for the whole test. MAD is
// far less sensitive to exactly that failure mode.
function bdsStatistic(x, m, epsilonMultiplier) {
  var scale = robustScale(x);
  if (scale === 0 || isNaN(scale)) return NaN;
  var epsilon = epsilonMultiplier * scale;
  var C1 = correlationIntegral(x, 1, epsilon);
  var Cm = correlationIntegral(x, m, epsilon);
  if (isNaN(C1) || isNaN(Cm)) return NaN;
  return Cm - Math.pow(C1, m);
}

// =================================================================
// PEAKS-OVER-THRESHOLD: Generalized Pareto Distribution fit +
// exceedance clustering test
// -----------------------------------------------------------------
// Method of Moments (not maximum likelihood) - a deliberate choice,
// not a shortcut: Nuuk's real exceedance counts at a 90th-percentile
// threshold are typically only ~15-20 points, too few for a numerical
// MLE optimizer to reliably converge in this already-fragile client-
// side sandbox. MOM has a simple closed-form formula (Hosking & Wallis
// 1987) - always computable, no iterative failure mode, and is
// specifically documented as competitive with or better than MLE at
// small sample sizes for the GPD.
// =================================================================
function gpdMomentFit(exceedances) {
  var n = exceedances.length;
  if (n < 5) return {xi: NaN, sigma: NaN, mean: NaN, n: n};
  var mean = meanArr(exceedances);
  var variance = exceedances.reduce(function(s, v) { return s + (v - mean) * (v - mean); }, 0) / n;
  if (variance <= 0) return {xi: NaN, sigma: NaN, mean: mean, n: n};
  var xi = 0.5 * (1 - mean * mean / variance);
  var sigma = 0.5 * mean * (mean * mean / variance + 1);
  return {xi: xi, sigma: sigma, mean: mean, n: n};
}

// Runs declustering: standard POT preprocessing step. GPD theory
// implicitly assumes exceedances are roughly independent events - a
// sustained real event (a plume that doesn't clear for 2-3 months)
// otherwise gets counted as 2-3 separate "independent" observations,
// inflating the effective sample size with duplicated information.
// Merges exceedances separated by a gap of <= runLength non-exceeding
// months into a single cluster, keeping only each cluster's MAXIMUM
// excess as its one representative value. runLength=1 (adjacent
// months only) matches the same adjacency definition the clustering
// test already uses, for consistency.
function declusterExceedances(x, threshold, runLength) {
  runLength = runLength || 1;
  var n = x.length;
  var clusters = [];
  var currentCluster = null;
  var lastExceedIdx = null;

  for (var i = 0; i < n; i++) {
    if (x[i] > threshold) {
      if (currentCluster !== null && (i - lastExceedIdx) <= runLength) {
        currentCluster.push(x[i] - threshold);
      } else {
        if (currentCluster !== null) clusters.push(currentCluster);
        currentCluster = [x[i] - threshold];
      }
      lastExceedIdx = i;
    }
  }
  if (currentCluster !== null) clusters.push(currentCluster);

  var clusterMaxima = clusters.map(function(c) { return Math.max.apply(null, c); });
  var clusterSizes = clusters.map(function(c) { return c.length; });
  return {clusterMaxima: clusterMaxima, nClusters: clusters.length, clusterSizes: clusterSizes};
}

// IMPORTANT LIMIT, stated honestly: IAAFT surrogates preserve the
// exact VALUE distribution of the real data (same sorted values, just
// reordered) - so the total COUNT of exceedances above any fixed
// threshold is mathematically GUARANTEED identical between the real
// series and every surrogate. That is not a testable quantity, and
// this tool does not pretend otherwise. What genuinely IS testable is
// CLUSTERING - whether extreme months land adjacent to each other in
// time more often than a random reordering of the same values would
// produce. This is the one part of extreme-value behavior that
// actually depends on temporal order, which is exactly what IAAFT
// surrogates scramble while holding everything else fixed.
function countExceedanceClusters(x, threshold) {
  var n = x.length;
  var pairCount = 0, totalExceed = 0;
  for (var i = 0; i < n; i++) {
    if (x[i] > threshold) {
      totalExceed++;
      if (i > 0 && x[i - 1] > threshold) pairCount++;
    }
  }
  return {pairCount: pairCount, totalExceed: totalExceed};
}

// =================================================================
// RECENT-vs-HISTORICAL comparison tests
// -----------------------------------------------------------------
// Both tests below ask a question about WHERE in time something
// happens - the same class of question as exceedance clustering
// above. Both use a FULL RANDOM PERMUTATION as the null, never IAAFT,
// for the same reason established with clustering: IAAFT's null
// hypothesis preserves exactly the kind of positional/low-frequency
// structure a real recent-vs-historical difference would show up as,
// biasing the test toward "not significant" even for a real shift.
// =================================================================

// Permutation test for a real difference between the early half and
// late half of an already-computed window-level trajectory (e.g. the
// AR(1) lambda trajectory) - reuses the trend test's own windowed
// estimates rather than re-fitting anything on a thin, unstable
// subset of the raw monthly data.
function permutationEarlyLateTest(validVals, nPerm, seed) {
  var n = validVals.length;
  if (n < 4) return {error: 'Only ' + n + ' valid windows - need at least 4 to compare early vs late.'};
  var half = Math.floor(n / 2);
  function meanDiff(arr) { return meanArr(arr.slice(half)) - meanArr(arr.slice(0, half)); }
  var realDiff = meanDiff(validVals);

  var rng = mulberry32(seed);
  var extreme = 0;
  var permDiffs = [];
  for (var p = 0; p < nPerm; p++) {
    var shuffled = validVals.slice();
    for (var i = shuffled.length - 1; i > 0; i--) {
      var j = Math.floor(rng() * (i + 1));
      var t = shuffled[i]; shuffled[i] = shuffled[j]; shuffled[j] = t;
    }
    var d = meanDiff(shuffled);
    permDiffs.push(d);
    if (Math.abs(d) >= Math.abs(realDiff)) extreme++;
  }
  return {realDiff: realDiff, pValue: extreme / nPerm, nEarly: half, nLate: n - half, permDiffs: permDiffs};
}

// Permutation test for whether exceedances are disproportionately
// concentrated in a recent window vs scattered evenly across the
// record. One-sided by design (>=, not |diff|) - the natural question
// is "did it get WORSE recently," not "did it change either way."
function testRecencyConcentration(exceedanceIndicator, recentStartIdx, nPerm, seed) {
  var n = exceedanceIndicator.length;
  var realRecentCount = 0;
  for (var i = recentStartIdx; i < n; i++) if (exceedanceIndicator[i]) realRecentCount++;

  var rng = mulberry32(seed);
  var extreme = 0;
  var permCounts = [];
  for (var p = 0; p < nPerm; p++) {
    var shuffled = exceedanceIndicator.slice();
    for (var i = shuffled.length - 1; i > 0; i--) {
      var j = Math.floor(rng() * (i + 1));
      var t = shuffled[i]; shuffled[i] = shuffled[j]; shuffled[j] = t;
    }
    var permRecentCount = 0;
    for (var i = recentStartIdx; i < n; i++) if (shuffled[i]) permRecentCount++;
    permCounts.push(permRecentCount);
    if (permRecentCount >= realRecentCount) extreme++;
  }
  return {realRecentCount: realRecentCount, recentMonths: n - recentStartIdx, pValue: extreme / nPerm, permCounts: permCounts};
}

// =================================================================
// EVENT COINCIDENCE ANALYSIS (cross-variable)
// -----------------------------------------------------------------
// Adapted from Donges, Schleussner, Siegmund & Donner (2016), "Event
// coincidence analysis for quantifying statistical interrelationships
// between event time series," Eur. Phys. J. Spec. Top. 225, 471-487 -
// a real, published method, applied in the literature to things like
// flood events preceding epidemic outbreaks and ENSO phases coinciding
// with precipitation extremes. Measures how often extreme months in
// ONE variable coincide (within a small time tolerance) with extreme
// months in ANOTHER variable.
//
// One deliberate departure from the original paper: Donges et al. get
// significance from an analytical formula assuming Poisson-process
// event statistics. This tool uses full random permutation instead,
// for the same small-sample reasoning applied everywhere else here
// (BDS, clustering, recency) - Nuuk's real exceedance counts (14-20)
// are far outside where an asymptotic Poisson approximation can be
// trusted, and permutation makes no such assumption.
// =================================================================

// Real coincidence rate: of variable A's extreme months, what
// fraction have an extreme month in variable B within +/- tau months?
function computeCoincidenceRate(indicatorA, indicatorB, tau) {
  var n = indicatorA.length;
  var aExtreme = [];
  for (var i = 0; i < n; i++) if (indicatorA[i]) aExtreme.push(i);
  if (aExtreme.length === 0) return {rate: NaN, nA: 0, nCoincide: 0};

  var nCoincide = 0;
  for (var k = 0; k < aExtreme.length; k++) {
    var idx = aExtreme[k];
    var found = false;
    for (var offset = -tau; offset <= tau; offset++) {
      var j = idx + offset;
      if (j >= 0 && j < n && indicatorB[j]) { found = true; break; }
    }
    if (found) nCoincide++;
  }
  return {rate: nCoincide / aExtreme.length, nA: aExtreme.length, nCoincide: nCoincide};
}

// Permutation test: reshuffle variable B's extreme-month positions
// (same total count, randomized WHERE they land), recompute the
// coincidence rate against the REAL, fixed variable A each time.
// Tests "would this much coincidence happen if B's extremes were
// scattered independently of A's actual extreme months."
function eventCoincidenceTest(indicatorA, indicatorB, tau, nPerm, seed) {
  var real = computeCoincidenceRate(indicatorA, indicatorB, tau);
  if (isNaN(real.rate)) return {error: 'No extreme months found in the reference variable.'};

  var rng = mulberry32(seed);
  var extreme = 0;
  var permRates = [];
  for (var p = 0; p < nPerm; p++) {
    var shuffledB = indicatorB.slice();
    for (var i = shuffledB.length - 1; i > 0; i--) {
      var j = Math.floor(rng() * (i + 1));
      var t = shuffledB[i]; shuffledB[i] = shuffledB[j]; shuffledB[j] = t;
    }
    var permResult = computeCoincidenceRate(indicatorA, shuffledB, tau);
    if (!isNaN(permResult.rate)) {
      permRates.push(permResult.rate);
      if (permResult.rate >= real.rate) extreme++;
    }
  }
  var pValue = permRates.length > 0 ? extreme / permRates.length : NaN;
  return {realRate: real.rate, nA: real.nA, nCoincide: real.nCoincide, pValue: pValue, permRates: permRates};
}

function corrAtLagFromAC1(ac1a, ac1b, lag) {
  // Takes ALREADY-COMPUTED AC1 trajectories and just lag-shifts + correlates.
  // PERF FIX: the previous version recomputed rollingAC1Arr() from scratch
  // inside every one of the 13 lag checks, for every one of 300 surrogates
  // (300 x 13 x 2 = 7,800 redundant AC1 passes instead of 300 x 2 = 600).
  // That ~13x waste, combined with a browser JS engine slower than local
  // Node, is the most likely cause of multi-minute runs. AC1 doesn't
  // depend on lag - only the shift-and-correlate step does - so it now
  // computes AC1 exactly once per surrogate.
  var x, y;
  if (lag === 0) { x = ac1a; y = ac1b; }
  else if (lag > 0) { x = ac1a.slice(0, ac1a.length - lag); y = ac1b.slice(lag); }
  else { var L = -lag; x = ac1a.slice(L); y = ac1b.slice(0, ac1b.length - L); }
  var validX = [], validY = [];
  for (var i=0;i<x.length;i++) {
    if (!isNaN(x[i]) && !isNaN(y[i])) { validX.push(x[i]); validY.push(y[i]); }
  }
  if (validX.length < 4) return NaN;
  return pearsonCorr(validX, validY);
}

function lagScanMaxStat(a, b, lags, window) {
  var ac1a = rollingAC1Arr(a, window);
  var ac1b = rollingAC1Arr(b, window);
  var corrs = {}, best = null, bestAbs = -1;
  lags.forEach(function(lag) {
    var c = corrAtLagFromAC1(ac1a, ac1b, lag);
    corrs[lag] = c;
    if (!isNaN(c) && Math.abs(c) > bestAbs) { bestAbs = Math.abs(c); best = lag; }
  });
  return {corrs: corrs, bestLag: best, maxAbs: bestAbs < 0 ? NaN : bestAbs};
}

// Full 4-step pipeline: preprocess -> lag-scan -> IAAFT surrogates -> p-value
// FIX: real browser testing confirmed the surrogate loop below is
// genuinely much slower in GEE's sandboxed JS environment than local
// benchmarks suggested (Node: ~650ms for 100 surrogates; real browser:
// slow enough to trigger Chrome's "Page Unresponsive" watchdog, which
// fires after several seconds of continuously blocked main thread).
// Running all 100 surrogates in one uninterrupted for-loop was the
// actual bottleneck. Rewritten as a CHUNKED, async, callback-based
// version: processes a small batch (default 5) of surrogates, then
// yields via ui.util.setTimeout before continuing with the next batch
// - the browser gets a real chance to paint/respond between batches,
// so no single tick blocks long enough to look "unresponsive," even
// if the total wall-clock time is similar or slightly longer overall.
function fullCouplingTestChunked(datesA, valsA, datesB, valsB, opts, onProgress, onComplete) {
  opts = opts || {};
  var window = opts.window || 24;
  var nSurr = opts.nSurr || 100;
  var lags = opts.lags || [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6];
  var seed = opts.seed || Math.floor(Math.random()*1e6);
  var chunkSize = opts.chunkSize || 1;

  var aligned = alignSeries(datesA, valsA, datesB, valsB);
  var n = aligned.dates.length;
  if (n < window + 12) {
    onComplete({error: 'Only ' + n + ' overlapping real months - need at least ' + (window+12) + ' for a window=' + window + ' AC1 test.'});
    return;
  }
  var aDes = deseasonalize(aligned.dates, aligned.a);
  var bDes = deseasonalize(aligned.dates, aligned.b);
  var aDet = linearDetrend(aDes);
  var bDet = linearDetrend(bDes);

  var real = lagScanMaxStat(aDet, bDet, lags, window);
  if (isNaN(real.maxAbs)) {
    onComplete({error: 'Could not compute a valid lag-scan statistic (too few AC1 windows for window=' + window + ').'});
    return;
  }

  var extreme = 0, validSurr = 0, surrStats = [];
  var s = 0;

  function processChunk() {
    var chunkStart = Date.now();
    var end = Math.min(s + chunkSize, nSurr);
    for (; s < end; s++) {
      var aSurr = iaaftSurrogate(aDet, 20, seed*10000+s);
      var bSurr = iaaftSurrogate(bDet, 20, seed*10000+s+50000);
      var surrRes = lagScanMaxStat(aSurr, bSurr, lags, window);
      if (isNaN(surrRes.maxAbs)) continue;
      validSurr++;
      surrStats.push(surrRes.maxAbs);
      if (surrRes.maxAbs >= real.maxAbs) extreme++;
    }
    var computeMs = Date.now() - chunkStart;
    if (onProgress) onProgress(s, nSurr);
    if (s < nSurr) {
      // DIAGNOSTIC: measures the gap between scheduling this callback and
      // it actually firing, isolating ui.util.setTimeout's OWN overhead
      // from the actual surrogate computation above. If "scheduling
      // overhead" below is large (many seconds) while "compute" stays
      // small, that confirms setTimeout itself - not the math, not
      // widget updates - is the real cost, and argues for fewer, larger
      // chunks rather than more, smaller ones.
      var scheduledAt = Date.now();
      ui.util.setTimeout(function() {
        var overheadMs = Date.now() - scheduledAt;
        print('[coupling] chunk done: compute=' + computeMs + 'ms, setTimeout scheduling overhead=' + overheadMs + 'ms');
        processChunk();
      }, 0);
    } else {
      var pValue = validSurr > 0 ? extreme / validSurr : NaN;
      onComplete({
        n: n, window: window, nSurr: validSurr, lags: lags,
        corrs: real.corrs, bestLag: real.bestLag, maxAbs: real.maxAbs,
        pValue: pValue, surrStats: surrStats
      });
    }
  }
  processChunk();
}

// =================================================================
// Chunked pipeline for the AR(1) relaxation-rate + BDS residual check.
// For EACH surrogate: generate one IAAFT fake series, then derive BOTH
// test statistics from that SAME surrogate (its own lambda-trend tau,
// AND its own AR(1) residuals' BDS stat) - reuses the expensive
// surrogate-generation step across both tests instead of drawing two
// independent surrogate sets, which is both cheaper and more
// statistically coherent (both tests are evaluated against the same
// realization of the null hypothesis, not two different ones).
// Same safe chunking pattern as fullCouplingTestChunked: chunkSize=1,
// yields via ui.util.setTimeout every single surrogate, minimal widget
// mutation during the loop (frequent print() instead).
// =================================================================
function modelFitAnalysisChunked(dates, vals, opts, onProgress, onComplete) {
  opts = opts || {};
  var window = opts.window || 36; // validated via direct simulation: window=24 (borrowed from the simpler AC1 test) was swamped by per-window phi estimation noise (std~0.18 around a true value of 0.5, at n=24). window=36 reliably detects real trends (18/20 test trials) while leaving enough windows for the trend test itself.
  var nSurr = opts.nSurr || 20;
  var nIter = opts.nIter || 20;
  var bdsM = opts.bdsM || 2;
  var bdsEpsMult = opts.bdsEpsMult || 0.7;
  var seed = opts.seed || Math.floor(Math.random() * 1e6);
  var chunkSize = 1;

  if (dates.length < window + 12) {
    onComplete({error: 'Only ' + dates.length + ' real months - need at least ' + (window + 12) + ' for a window=' + window + ' model fit.'});
    return;
  }

  var des = deseasonalize(dates, vals);
  var det = linearDetrend(des);

  var realLambdas = rollingLambdaTrajectory(det, window);
  var validLambdas = realLambdas.filter(function(v) { return !isNaN(v); });
  if (validLambdas.length < 4) {
    onComplete({error: 'Only ' + validLambdas.length + ' windows produced a valid (stationary, phi in 0..1) relaxation rate - too few to test a trend. Try a shorter window.'});
    return;
  }
  var realTau = kendallTau(validLambdas);

  var realFit = fitAR1(det);
  var robustFit = fitAR1Robust(det);
  if (isNaN(realFit.phi)) {
    onComplete({error: 'Could not fit a valid AR(1) model to the whole series.'});
    return;
  }
  var realBDS = bdsStatistic(realFit.residuals, bdsM, bdsEpsMult);
  if (isNaN(realBDS)) {
    onComplete({error: 'Could not compute a valid BDS statistic on the residuals (series may be too short or too uniform).'});
    return;
  }

  var s = 0;
  var tauExtreme = 0, tauValidSurr = 0, tauSurrStats = [];
  var bdsExtreme = 0, bdsValidSurr = 0, bdsSurrStats = [];

  function processChunk() {
    var chunkStart = Date.now();
    var end = Math.min(s + chunkSize, nSurr);
    for (; s < end; s++) {
      var surr = iaaftSurrogate(det, nIter, seed * 10000 + s);

      var surrLambdas = rollingLambdaTrajectory(surr, window);
      var surrValidLambdas = surrLambdas.filter(function(v) { return !isNaN(v); });
      if (surrValidLambdas.length >= 4) {
        var surrTau = kendallTau(surrValidLambdas);
        if (!isNaN(surrTau)) {
          tauValidSurr++;
          tauSurrStats.push(surrTau);
          if (Math.abs(surrTau) >= Math.abs(realTau)) tauExtreme++;
        }
      }

      var surrFit = fitAR1(surr);
      if (!isNaN(surrFit.phi)) {
        var surrBDS = bdsStatistic(surrFit.residuals, bdsM, bdsEpsMult);
        if (!isNaN(surrBDS)) {
          bdsValidSurr++;
          bdsSurrStats.push(surrBDS);
          if (Math.abs(surrBDS) >= Math.abs(realBDS)) bdsExtreme++;
        }
      }
    }
    var computeMs = Date.now() - chunkStart;
    if (onProgress) onProgress(s, nSurr);
    if (s < nSurr) {
      ui.util.setTimeout(function() {
        if (s % 10 === 0) print('[model-fit] progress: ' + s + '/' + nSurr + ' (last chunk compute=' + computeMs + 'ms)');
        processChunk();
      }, 0);
    } else {
      var tauPValue = tauValidSurr > 0 ? tauExtreme / tauValidSurr : NaN;
      var bdsPValue = bdsValidSurr > 0 ? bdsExtreme / bdsValidSurr : NaN;
      onComplete({
        n: dates.length, window: window, nSurr: nSurr,
        realLambdas: realLambdas, realTau: realTau, tauPValue: tauPValue, tauNSurr: tauValidSurr, tauSurrStats: tauSurrStats,
        phi: realFit.phi, phiRobust: robustFit.phi, residuals: realFit.residuals, realBDS: realBDS, bdsPValue: bdsPValue, bdsNSurr: bdsValidSurr, bdsSurrStats: bdsSurrStats
      });
    }
  }
  processChunk();
}

// =================================================================
// Chunked pipeline for peaks-over-threshold analysis. The GPD fit
// itself is NOT tested against surrogates (it can't be - IAAFT
// surrogates share the real data's exact value distribution by
// construction, so GPD parameters would come out statistically
// identical regardless). What IS genuinely tested: whether the real
// exceedances cluster together in time more than chance reordering
// (via IAAFT) of the same values would produce.
// =================================================================
function extremeValueAnalysisChunked(dates, vals, opts, onProgress, onComplete) {
  opts = opts || {};
  var percentile = opts.percentile || 90;
  var nSurr = opts.nSurr || 20;
  var seed = opts.seed || Math.floor(Math.random() * 1e6);
  var chunkSize = 1;

  if (dates.length < 30) {
    onComplete({error: 'Only ' + dates.length + ' real months - need at least 30 for a meaningful exceedance analysis.'});
    return;
  }

  var des = deseasonalize(dates, vals);
  var det = linearDetrend(des);

  var sorted = det.slice().sort(function(a, b) { return a - b; });
  var idx = Math.min(Math.floor((percentile / 100) * sorted.length), sorted.length - 1);
  var threshold = sorted[idx];

  var realExceedances = [];
  det.forEach(function(v) { if (v > threshold) realExceedances.push(v - threshold); });
  if (realExceedances.length < 5) {
    onComplete({error: 'Only ' + realExceedances.length + ' exceedances above the ' + percentile + 'th percentile - too few for a reliable GPD fit. Try a lower percentile.'});
    return;
  }

  var realGPD = gpdMomentFit(realExceedances);
  var realClusters = countExceedanceClusters(det, threshold);

  // Declustered fit: merges adjacent exceedances into single real
  // events first (runs declustering, run length 1 month, matching the
  // clustering test's own adjacency definition), then refits GPD to
  // just the cluster maxima. Reported alongside the raw fit, not in
  // place of it, so the difference itself is visible.
  var declusterResult = declusterExceedances(det, threshold, 1);
  var declusteredGPD = declusterResult.clusterMaxima.length >= 5 ? gpdMomentFit(declusterResult.clusterMaxima) : null;

  // Recency-concentration setup: "recent" = actual last 5 calendar
  // years, computed from real dates (not array-index fraction, which
  // would be skewed by real gaps in the record). Manual loop, not
  // Array.prototype.findIndex - this sandbox is already known to lack
  // other ES6 array methods (Array.fill), so avoid relying on it here.
  var lastDate = dates[dates.length - 1];
  var lastYear = parseInt(lastDate.slice(0, 4), 10);
  var lastMonth = lastDate.slice(5, 7);
  var cutoffDateStr = (lastYear - 5) + '-' + lastMonth;
  var recentStartIdx = -1;
  for (var di = 0; di < dates.length; di++) {
    if (dates[di] >= cutoffDateStr) { recentStartIdx = di; break; }
  }
  if (recentStartIdx === -1) recentStartIdx = Math.floor(dates.length / 2);
  var exceedanceIndicator = det.map(function(v) { return v > threshold ? 1 : 0; });
  var recencyResult = (recentStartIdx > 0 && dates.length - recentStartIdx >= 6) ?
    testRecencyConcentration(exceedanceIndicator, recentStartIdx, 1000, seed + 555) : null;

  // FIX (caught by direct testing, not assumed): the clustering null
  // hypothesis below uses a FULL RANDOM PERMUTATION, not an IAAFT
  // surrogate - a real, substantive difference, not just a style
  // choice. A direct test with an obvious 8-month cluster injected
  // showed IAAFT surrogates gave pairCounts of 12-17 (HIGHER than the
  // real cluster's 8!) because IAAFT's null hypothesis specifically
  // preserves the same autocorrelation/spectral structure that real
  // clustering creates - it is the wrong null for testing clustering
  // itself. A full permutation destroys ALL temporal order, which is
  // the standard, correct null hypothesis for extremal clustering
  // tests in the literature. Confirmed by the same test: permutation
  // correctly gave pairCounts of 0-4, p=0.000, for the same injected
  // cluster where IAAFT gave p=1.000 (backwards).
  function fisherYatesShuffle(x, drawSeed) {
    var rng = mulberry32(drawSeed);
    var a = x.slice();
    for (var i = a.length - 1; i > 0; i--) {
      var j = Math.floor(rng() * (i + 1));
      var t = a[i]; a[i] = a[j]; a[j] = t;
    }
    return a;
  }

  var s = 0;
  var clusterExtreme = 0, clusterValidSurr = 0, clusterSurrStats = [];

  function processChunk() {
    var chunkStart = Date.now();
    var end = Math.min(s + chunkSize, nSurr);
    for (; s < end; s++) {
      var shuffled = fisherYatesShuffle(det, seed * 10000 + s);
      var surrClusters = countExceedanceClusters(shuffled, threshold);
      clusterValidSurr++;
      clusterSurrStats.push(surrClusters.pairCount);
      if (surrClusters.pairCount >= realClusters.pairCount) clusterExtreme++;
    }
    var computeMs = Date.now() - chunkStart;
    if (onProgress) onProgress(s, nSurr);
    if (s < nSurr) {
      ui.util.setTimeout(function() {
        if (s % 10 === 0) print('[extreme-value] progress: ' + s + '/' + nSurr + ' (compute=' + computeMs + 'ms)');
        processChunk();
      }, 0);
    } else {
      var clusterPValue = clusterValidSurr > 0 ? clusterExtreme / clusterValidSurr : NaN;
      onComplete({
        n: dates.length, percentile: percentile, threshold: threshold,
        nExceedances: realExceedances.length, gpdXi: realGPD.xi, gpdSigma: realGPD.sigma, meanExcess: realGPD.mean,
        nClusters: declusterResult.nClusters, declGpdXi: declusteredGPD ? declusteredGPD.xi : NaN, declGpdSigma: declusteredGPD ? declusteredGPD.sigma : NaN, declGpdMean: declusteredGPD ? declusteredGPD.mean : NaN,
        realPairCount: realClusters.pairCount, totalExceed: realClusters.totalExceed,
        clusterPValue: clusterPValue, clusterNSurr: clusterValidSurr, clusterSurrStats: clusterSurrStats,
        recency: recencyResult
      });
    }
  }
  processChunk();
}

// surface = <20 dbar, mid = 20-300 dbar, deep = >300 dbar. Monthly means, real GEM MarineBasis CTD casts.
var NUUK_DATA = {
  turbidity: {
    surface: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-07','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.7808,0.7746,0.5834,0.6936,1.1939,0.7198,0.5878,0.7719,0.5112,0.6236,0.9653,0.6311,0.6903,0.5768,0.6612,0.8559,0.8763,0.9522,0.6746,0.688,0.7872,0.8433,0.7284,0.6244,0.8139,0.6952,0.5638,0.6756,0.446,0.9931,0.7783,0.9147,0.8832,0.8838,0.832,0.6719,0.6082,1.0012,0.7676,0.7001,0.712,0.77,0.6841,0.9289,0.7746,0.8466,0.6477,0.6246,1.0299,1.6621,1.706,1.7182,0.8352,0.8064,0.7166,0.6042,0.8117,0.8941,1.044,0.8106,0.7851,0.8115,0.6652,0.8315,0.7348,0.5377,1.0462,0.8828,0.7739,0.8717,0.6991,1.1167,1.0157,0.7482,0.8217,0.6764,0.6788,0.526,0.6673,0.7112,0.5694,0.4767,0.7728,0.7533,0.9823,0.6509,0.6967,0.6833,0.5483,0.9041,0.6952,0.7515,0.6316,0.9359,0.6729,0.6362,0.7213,0.5755,0.6091,0.5623,0.4889,0.615,0.5469,0.5929,0.5602,0.5765,0.5912,0.4974,0.7344,0.5149,0.7035,0.7481,0.8748,0.4446,0.7146,0.4895,0.5505,0.5422,0.3721,0.4962,0.639,0.5484,0.5443,0.5636,1.1517,1.4734,1.0865,2.1216,0.9026,1.0819,0.6528,0.5327,0.4597,0.5377,0.5454,0.708,0.5117,0.4449,0.5838,0.4148,0.5341,0.8671,1.1179,0.5625,0.4909,0.5984,0.4956,0.7658,0.4829,0.3799,0.5691,0.4549,0.5369,0.8137,0.5049,0.4922,0.5552,0.5772,0.4891,0.5175,0.4772,0.3501,0.6092,0.7692,0.6513,0.4161,0.8087,0.6254,0.7875,0.3935,0.7706,0.5694,1.3637,2.39,3.2433,0.4737,0.3288,0.3889,0.2266,0.5373,0.6914,0.4752,0.3589,0.4124,0.9951,1.3779,1.5952,0.0827,2.9577,0.1443,0.4136,0.3105,0.3237,0.3122,0.2055]},
    mid: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-07','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.5951,0.5123,0.4762,0.593,0.4957,0.6804,0.5975,0.612,0.4934,0.5453,0.6538,0.648,0.6174,0.4754,0.7462,0.865,0.5499,0.6764,1.0392,1.1469,0.9091,0.7885,1.1062,0.4925,0.8546,1.1044,0.782,0.8983,0.6702,1.1959,0.7151,0.6707,0.9219,1.0089,0.6588,0.4868,0.53,0.6315,0.5015,0.4954,0.4541,0.6204,0.7182,0.7269,0.7806,0.7733,0.6503,0.6052,0.6624,1.5654,1.5621,1.6517,1.2601,0.6423,0.7293,0.5977,0.4873,0.4711,0.4907,0.5519,0.7575,0.8459,0.6269,0.5178,0.6697,0.5176,0.5336,0.6345,0.6463,0.9374,0.8219,0.672,0.7667,0.5318,0.6445,0.7217,0.5513,0.5623,0.4474,0.6769,0.5992,0.5343,0.7715,0.7767,0.7414,0.6874,0.6309,0.604,0.5705,0.486,0.4937,0.5323,0.7354,0.6194,0.6152,0.5461,0.5076,0.5803,0.6318,0.4366,0.3929,0.5786,0.5858,0.6377,0.461,0.6781,0.5431,0.5065,0.5892,0.5286,0.5915,0.5357,0.4919,0.4888,0.4804,0.5341,0.5478,0.4549,0.4066,0.3976,0.3772,0.4375,0.5094,0.6155,0.8712,1.1021,1.1525,1.9001,0.9723,0.9372,0.3929,0.3806,0.4737,0.391,0.5684,0.5352,0.462,0.4777,0.5598,0.4963,0.4236,0.6074,0.9881,0.443,0.5018,0.44,0.5482,0.4543,0.4622,0.3878,0.3554,0.4215,0.5195,0.4129,0.4298,0.6429,0.5142,0.5747,0.3382,0.4641,0.4473,0.3379,0.4873,0.455,0.281,0.4437,0.4825,0.6795,0.4425,0.4424,0.3734,0.6656,0.6804,1.4118,2.8208,0.524,0.3009,0.4059,0.3423,0.3913,0.2844,0.3608,0.4229,0.5254,0.486,0.5349,1.674,0.9303,0.9649,0.3383,0.435,0.2433,0.3757,0.295,0.234]},
    deep: {dates: ['2005-10','2005-11','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2008-02','2008-03','2008-05','2008-06','2008-07','2008-08','2009-01','2009-04','2009-05','2009-06','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-07','2010-12','2011-01','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-12','2012-01','2012-04','2012-05','2012-06','2012-07','2012-09','2012-10','2012-12','2013-01','2013-03','2013-04','2013-05','2013-06','2013-09','2013-10','2013-11','2013-12','2014-01','2014-03','2014-04','2014-06','2014-07','2014-08','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-08','2019-09','2019-11','2020-01','2020-05','2020-06','2020-07','2020-08','2020-11','2021-01','2021-02','2021-03','2021-05','2021-06','2021-07','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-11','2023-12','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.7173,0.528,0.6643,0.5038,0.5576,0.7515,0.853,0.6372,0.4394,0.8999,0.9518,0.5002,0.6906,1.4303,1.4636,1.3411,1.2092,1.14,0.9951,0.9468,1.3024,0.847,1.5317,0.563,0.5935,0.4586,0.6776,0.7587,0.6992,0.7874,0.7534,0.6493,0.5869,0.7065,1.5387,1.5122,2.3545,0.599,0.8029,0.7306,0.6029,0.4799,0.7429,0.838,0.473,0.6884,0.5963,0.5454,0.8974,1.3595,1.03,0.6413,0.526,0.6489,0.7513,0.5448,0.6697,0.6407,0.5599,0.9449,0.8245,0.7185,0.6245,0.5776,0.5437,0.651,0.5717,0.4887,0.5333,0.5702,0.7,0.4748,0.4427,0.8467,0.6142,0.7093,0.742,0.5522,0.523,0.633,0.5799,0.6982,0.7146,0.487,0.5584,0.569,0.4759,0.4715,0.4601,0.4383,0.4745,0.6725,0.8978,1.0908,1.1674,1.913,1.0925,0.9631,0.4587,0.5859,0.3941,0.5683,0.5562,0.4436,0.494,0.5317,0.6403,0.8765,0.7341,0.5704,0.6511,0.5667,0.5185,0.4053,0.5236,0.1371,0.5242,0.4504,0.3118,0.4655,0.361,0.5449,0.5505,0.4738,0.6952,0.4025,0.5028,0.3411,0.4792,0.8141,0.2814,0.4729,0.3045,0.4181,0.4843,0.5134,0.3821,0.6814,0.4793,0.1413,0.0335,0.508,0.6108,0.5344,0.6033,0.6338,0.4226]},
    midUpper: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-07','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.7138,0.5955,0.4998,0.6399,0.7608,0.6745,0.5527,0.6978,0.5061,0.5512,0.7523,0.584,0.6226,0.5501,0.6647,0.8681,0.6685,0.781,0.7025,0.5932,0.6016,0.7293,0.7455,0.5607,0.6414,0.693,0.7069,0.8144,0.4512,0.8434,0.6704,0.6765,0.7987,0.9073,0.6771,0.5346,0.516,0.8647,0.6832,0.5156,0.5123,0.5993,0.6492,0.8263,0.7606,0.7882,0.5419,0.6078,0.7754,1.5976,1.6247,1.6586,0.8051,0.7516,0.6949,0.5729,0.6756,0.6144,0.495,0.6051,0.7381,0.8291,0.5793,0.5711,0.6895,0.4934,0.5592,0.7279,0.8027,0.8492,0.7094,0.8306,0.8561,0.5241,0.7489,0.7285,0.6154,0.5005,0.4961,0.6953,0.5688,0.5284,0.7418,0.7609,0.7593,0.6721,0.6412,0.6521,0.5554,0.5804,0.4739,0.5482,0.6627,0.5766,0.6332,0.6154,0.5863,0.5905,0.6248,0.4538,0.3856,0.6349,0.5505,0.5651,0.4029,0.6024,0.5845,0.4972,0.6111,0.4995,0.5661,0.5528,0.3838,0.4388,0.5637,0.5022,0.5581,0.4579,0.3547,0.4074,0.3985,0.4573,0.5104,0.5878,0.9855,1.2443,1.1065,1.9908,0.9235,0.9697,0.4266,0.4327,0.4446,0.3785,0.5222,0.5969,0.4569,0.4829,0.5734,0.4368,0.438,0.7124,1.0702,0.5176,0.436,0.4628,0.512,0.4572,0.4471,0.3993,0.4856,0.4224,0.371,0.361,0.3523,0.5071,0.5223,0.8841,0.3375,0.4048,0.4942,0.3219,0.5629,0.5775,0.3346,0.4063,0.6336,0.6537,0.5232,0.4019,0.4922,1.9915,1.5558,2.2919,3.367,0.5026,0.3068,0.3781,0.2305,0.3872,0.2927,0.3513,0.3292,0.4056,0.5726,0.799,2.5687,0.3565,3.0244,0.1992,0.3936,0.2438,0.2511,0.2757,0.1603]},
    midCore: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-07','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.6161,0.5023,0.475,0.6009,0.4494,0.6325,0.5389,0.5822,0.4958,0.5627,0.6201,0.5333,0.5799,0.4879,0.7111,0.8603,0.5966,0.6745,0.8723,0.9745,0.6306,0.7349,0.9111,0.4805,0.6937,1.3416,0.8082,0.8766,0.5787,1.0032,0.6446,0.6453,0.8754,0.976,0.6306,0.5229,0.5126,0.6379,0.5173,0.4768,0.426,0.5901,0.6752,0.7317,0.7819,0.7828,0.6529,0.6033,0.6254,1.5743,1.5714,1.6348,0.8066,0.6233,0.7069,0.5886,0.44,0.4596,0.4396,0.5069,0.7536,0.853,0.6332,0.5503,0.6786,0.4928,0.5005,0.6016,0.5947,0.8616,0.7494,0.7063,0.7526,0.532,0.6412,0.7115,0.5717,0.4935,0.4184,0.6946,0.5441,0.5656,0.7194,0.7614,0.7345,0.6593,0.6407,0.6173,0.5568,0.4536,0.4421,0.5139,0.6828,0.6463,0.6006,0.5836,0.5002,0.5925,0.6016,0.4223,0.3701,0.5482,0.5692,0.6126,0.4656,0.6544,0.5492,0.513,0.5643,0.498,0.5484,0.5055,0.4569,0.4298,0.4783,0.5381,0.5347,0.4559,0.3744,0.373,0.3578,0.4356,0.4833,0.5938,0.8444,1.0798,1.1418,1.8732,0.9447,0.9117,0.3783,0.3475,0.4493,0.3955,0.5587,0.5209,0.4636,0.4581,0.5491,0.4462,0.3813,0.5712,1.0085,0.4518,0.4766,0.3873,0.5238,0.4501,0.4449,0.3883,0.3458,0.3699,0.5557,0.3469,0.4383,0.6046,0.5034,0.7209,0.3585,0.5339,0.4798,0.3316,0.4615,0.4247,0.2986,0.4257,0.4709,0.6844,0.4508,0.4316,0.3765,0.5288,0.3571,1.4424,3.2505,0.5315,0.3226,0.3793,0.2934,0.3532,0.2318,0.3585,0.396,0.4864,0.5529,0.6373,1.9436,1.9822,1.7037,0.2666,0.3587,0.1836,0.2839,0.1968,0.1753]},
    midLower: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-07','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.5576,0.5025,0.4722,0.5783,0.4702,0.7151,0.6462,0.6148,0.4892,0.5328,0.6566,0.7366,0.6412,0.4523,0.7857,0.8676,0.4955,0.6569,1.3553,1.3712,1.1579,0.8692,1.3071,0.4805,1.0036,1.0291,0.7792,0.9292,0.7743,1.3935,0.8084,0.6866,0.9851,1.0508,0.6779,0.4511,0.5443,0.581,0.4549,0.5037,0.4616,0.6445,0.7603,0.704,0.7838,0.7641,0.6701,0.6059,0.6645,1.5499,1.5434,1.6633,1.6509,0.6332,0.7509,0.6087,0.4812,0.4502,0.525,0.5725,0.7639,0.8445,0.6323,0.4856,0.6596,0.5388,0.5504,0.6456,0.6495,1.005,0.9007,0.6177,0.7542,0.5333,0.626,0.7272,0.5251,0.639,0.4576,0.6562,0.6499,0.5147,0.812,0.79,0.7433,0.7091,0.6223,0.5851,0.5826,0.4886,0.5318,0.5488,0.7874,0.6079,0.6214,0.5075,0.4968,0.5703,0.6532,0.4427,0.4094,0.5874,0.6039,0.6687,0.4716,0.7088,0.5309,0.5041,0.6013,0.5545,0.6252,0.5554,0.5419,0.5388,0.4652,0.5378,0.5543,0.4536,0.4382,0.4119,0.3898,0.4348,0.5278,0.6354,0.8662,1.0886,1.1686,1.8998,1.0002,0.9476,0.3983,0.3921,0.5003,0.3905,0.5839,0.5324,0.4619,0.4896,0.5643,0.5413,0.4492,0.6105,0.9374,0.4149,0.5315,0.4783,0.5716,0.4578,0.4766,0.3852,0.3358,0.4555,0.5468,0.4704,0.4397,0.7154,0.5198,0.4052,0.3248,0.4229,0.416,0.3452,0.4895,0.4477,0.2516,0.4631,0.4602,0.6814,0.421,0.4575,0.3477,0.4927,0.7206,1.2167,2.322,0.5233,0.2854,0.429,0.3969,0.4172,0.3175,0.3895,0.5675,0.5751,0.4245,0.4061,0.6958,0.3438,0.0482,0.4135,0.4938,0.2828,0.4612,0.3638,0.2876]},
  },
  fluorescence: {
    surface: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [1.0603,0.7148,0.1478,0.0741,0.0547,0.0939,0.5051,3.5286,0.1722,0.8036,0.8088,0.5607,0.7944,0.1988,0.0826,0.1169,0.1658,4.8229,1.1436,1.1285,0.6823,1.1284,0.152,0.135,0.0762,0.1845,0.3561,0.7749,0.2809,0.5213,0.9149,1.2148,0.754,0.1396,0.0874,0.1046,0.204,0.1975,4.8993,0.3384,0.5801,1.5331,1.182,0.1792,0.2178,0.27,0.267,0.1649,0.1641,0.2364,0.7127,0.2509,1.6835,0.9428,0.7014,0.816,0.2034,0.4151,0.3915,0.4195,0.1765,2.7155,2.0194,0.4457,1.1053,1.17,0.7897,1.6791,0.1266,0.0723,0.0936,4.3008,1.6194,0.2018,1.5772,0.2725,0.7222,0.5273,0.1838,0.0984,0.1304,0.2609,3.7809,1.8029,1.4971,0.8676,0.9389,0.3713,0.3438,0.2317,0.1433,0.0996,0.0969,0.2004,1.1507,2.4211,0.2031,2.0847,1.5726,1.4367,0.2383,0.1446,0.0815,0.0729,0.0916,0.2939,2.8434,1.2452,0.8962,0.3124,0.3282,0.2544,0.1194,0.1329,0.1537,0.3481,1.0443,0.385,0.7942,1.2512,0.485,0.1175,0.0801,0.3193,0.1314,1.9234,0.2958,1.8205,0.6929,0.8425,0.2603,0.1104,0.0583,0.0976,0.6299,4.7143,0.4644,0.2675,1.1395,0.4815,0.2693,0.1476,0.0766,0.0818,0.2007,0.4058,2.0214,0.2835,0.9937,1.3659,1.2144,0.9283,0.2982,0.2239,0.0814,0.064,1.0995,0.1406,1.242,0.6961,1.4512,0.4396,0.2198,0.1512,0.1002,0.1166,0.3213,0.2308,1.0541,0.4861,0.4436,1.3694,0.6819,0.4773,0.1477,0.0768,0.0882,0.723,4.3236,0.9264,0.6858,0.3726,0.2328,0.1338,0.3507,0.9677,0.1376,0.722,1.0569,0.4942,0.7218,0.3,0.2678,0.1136,0.3528,1.6166,0.0753,0.0915,0.5419,0.5352,0.5137,0.2975,0.1982]},
    mid: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.2293,0.1528,0.1053,0.0612,0.0563,0.1384,0.4967,1.5125,0.2991,0.2978,0.1896,0.1451,0.1982,0.0944,0.0854,0.1724,0.2252,1.6208,0.7788,0.2613,0.1727,0.1618,0.1396,0.1225,0.07,0.1989,0.3671,0.9512,0.1151,0.1969,0.3522,0.2588,0.2155,0.0913,0.0853,0.1146,0.2277,0.3858,1.4613,0.1851,0.1491,0.1784,0.2064,0.1425,0.1093,0.255,0.2442,0.1552,0.2377,0.2367,0.6444,0.1726,0.3019,0.1505,0.2291,0.148,0.128,0.3728,0.3551,0.3769,0.3258,1.1217,1.8505,0.6372,0.2272,0.1438,0.0919,0.2256,0.1122,0.0773,0.1188,0.6191,1.3176,0.2531,0.1805,0.2368,0.1283,0.1782,0.1254,0.104,0.1252,0.4552,1.5903,2.1265,0.9098,0.3854,0.2689,0.1502,0.1488,0.1213,0.0928,0.1008,0.1214,0.2899,0.6264,2.1034,0.4586,0.5378,0.195,0.2593,0.1679,0.1166,0.0722,0.1199,0.1166,0.5311,3.0675,0.3865,0.217,0.2058,0.2501,0.1467,0.1015,0.1573,0.2382,0.5496,0.4415,0.4298,0.2665,0.1721,0.1453,0.1028,0.12,0.3495,0.4052,0.7438,0.1301,0.3241,0.1426,0.2616,0.1151,0.0923,0.0716,0.0912,0.7604,1.5322,0.61,0.1776,0.1556,0.1555,0.1742,0.1197,0.0882,0.1264,0.2704,0.4231,1.6885,0.2823,0.3974,0.2862,0.2623,0.2364,0.1498,0.1248,0.0886,0.0858,0.6199,0.2707,0.4721,0.252,0.3029,0.1526,0.1218,0.1218,0.0887,0.1131,0.3091,0.2339,0.8578,0.4507,0.2305,0.2748,0.0976,0.2185,0.1296,0.0734,0.1407,1.1989,0.6151,0.3821,0.2864,0.1282,0.1282,0.1655,0.3783,0.7566,0.2592,0.1844,0.1885,0.247,0.1741,0.1537,0.0979,0.1525,0.2908,0.3939,0.1006,0.0969,0.3214,0.1457,0.1858,0.14,0.1158]},
    deep: {dates: ['2005-10','2005-11','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2008-02','2008-03','2008-05','2008-06','2008-07','2008-08','2008-10','2009-01','2009-04','2009-05','2009-06','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-07','2010-08','2010-09','2010-12','2011-01','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-12','2012-01','2012-04','2012-05','2012-06','2012-07','2012-09','2012-10','2012-12','2013-01','2013-03','2013-04','2013-05','2013-06','2013-09','2013-10','2013-11','2013-12','2014-01','2014-03','2014-04','2014-06','2014-07','2014-08','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-08','2019-09','2019-11','2020-01','2020-05','2020-06','2020-07','2020-08','2020-11','2021-01','2021-02','2021-03','2021-05','2021-06','2021-07','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-11','2023-12','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.0753,0.0772,1.053,0.1865,0.1556,0.0867,0.0675,0.1058,0.0845,0.0912,0.1764,0.136,0.8872,0.2913,0.0889,0.0944,0.058,0.1617,0.7784,0.0785,0.0609,0.0825,0.1551,0.0801,0.5519,0.6894,0.1279,0.1131,0.0776,0.089,0.2467,0.2421,0.1718,0.2526,0.1864,0.3896,0.0811,0.0829,0.0813,0.3536,0.3535,0.31,1.1524,2.4442,0.1927,0.0569,0.0711,0.1144,0.0819,0.7061,0.5463,0.0949,0.0717,0.0639,0.0825,0.1034,0.0945,0.4375,0.9567,2.06,0.5325,0.0737,0.0657,0.0885,0.0851,0.1003,0.2424,0.5619,0.3282,0.0839,0.0747,0.0678,0.1254,0.0863,0.4831,3.4137,0.1563,0.1082,0.1197,0.1521,0.099,0.167,0.2551,0.3475,0.1913,0.1712,0.1372,0.0884,0.1546,0.3596,0.346,0.6287,0.0791,0.0704,0.0732,0.1102,0.0814,0.0539,0.0584,0.6448,1.1383,0.3013,0.1142,0.0956,0.1739,0.101,0.089,0.1433,0.2843,0.3042,1.1416,0.1526,0.133,0.1029,0.1178,0.0933,0.2806,0.1213,0.1225,0.0961,0.1122,0.0895,0.0784,0.2675,0.3257,0.2169,0.0935,0.0991,0.1099,0.0671,0.1225,1.0008,0.2992,0.1751,0.123,0.1031,0.1778,0.2466,0.2818,0.1119,0.0768,0.0837,0.1036,0.0591,0.2486,0.0984,0.0807,0.1442,0.0993,0.0824,0.1141,0.078]},
    midUpper: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.6176,0.4338,0.1355,0.0694,0.0627,0.1432,0.8378,3.209,0.2924,0.691,0.5826,0.3579,0.5702,0.1491,0.0831,0.1738,0.3204,3.5613,1.33,0.6931,0.4594,0.4201,0.1547,0.1349,0.0823,0.2348,0.4802,1.0613,0.177,0.5317,0.7981,0.811,0.4398,0.1124,0.0868,0.1372,0.4147,0.2725,5.1947,0.5036,0.4169,0.6696,0.4659,0.1944,0.1866,0.2672,0.2515,0.1539,0.2371,0.365,0.979,0.2183,1.2811,0.5501,0.6261,0.3865,0.1863,0.3943,0.3717,0.4149,0.3262,2.7452,2.0475,0.9153,1.0971,0.486,0.212,0.8805,0.1237,0.0758,0.1758,1.1273,2.2998,0.5467,0.6355,0.4878,0.4127,0.4701,0.1697,0.1108,0.1658,0.625,3.3825,2.215,2.2665,0.9881,0.9553,0.3322,0.3146,0.1918,0.1239,0.0998,0.1273,0.3835,1.1753,3.6676,0.4971,2.4645,0.7903,0.5649,0.2325,0.131,0.0754,0.1173,0.1778,0.5877,3.8166,0.8297,0.7732,0.424,0.4058,0.2438,0.1246,0.1489,0.2413,0.6161,0.9335,0.8113,0.8021,0.4168,0.2817,0.1071,0.0904,0.3446,0.3471,1.2458,0.2547,1.2946,0.3666,0.6651,0.1784,0.105,0.1002,0.1292,0.8138,3.7513,0.6664,0.3785,0.4111,0.3698,0.1965,0.1441,0.084,0.1033,0.2963,0.4797,2.8225,0.4505,1.6602,0.8426,0.6936,0.6399,0.2194,0.1702,0.0826,0.0893,1.224,0.346,2.5438,0.931,0.6524,0.2585,0.1495,0.1312,0.0941,0.1452,0.4354,0.4436,1.3761,0.7989,0.6596,0.7463,0.1864,0.4327,0.1614,0.0804,0.1956,1.18,2.4982,0.8577,0.5996,0.2349,0.1549,0.1698,0.6824,1.6182,0.2843,0.7053,0.7157,0.4938,0.4096,0.2794,0.1815,0.1558,0.461,0.701,0.0901,0.111,0.6433,0.3062,0.4416,0.2366,0.1814]},
    midCore: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.2903,0.1698,0.1141,0.0611,0.0542,0.1436,0.5839,1.6186,0.333,0.3403,0.2127,0.1932,0.2209,0.091,0.0829,0.1663,0.2592,1.7401,0.9529,0.2978,0.2053,0.1251,0.1441,0.1203,0.0782,0.2089,0.4111,1.0197,0.1224,0.2215,0.4324,0.2991,0.2249,0.0941,0.0871,0.1248,0.2155,0.2627,1.3983,0.181,0.1607,0.1601,0.2654,0.1703,0.1111,0.261,0.2446,0.1409,0.2354,0.2678,0.7403,0.174,0.3158,0.134,0.3005,0.1431,0.1575,0.3799,0.3527,0.3753,0.3377,0.7227,1.5767,0.8103,0.2012,0.1376,0.1004,0.2018,0.1072,0.0729,0.1445,0.3831,1.6082,0.3521,0.1934,0.2071,0.133,0.1954,0.1405,0.1077,0.1278,0.3829,1.7394,2.2316,1.0577,0.4425,0.2769,0.1696,0.1874,0.1469,0.0975,0.1023,0.1199,0.3131,0.5319,2.2322,0.4863,0.5284,0.1917,0.2736,0.1991,0.121,0.0747,0.1175,0.1271,0.5635,2.7747,0.5053,0.2001,0.2495,0.3219,0.1684,0.1075,0.1554,0.2329,0.6486,0.5109,0.6616,0.3114,0.184,0.1647,0.0993,0.0846,0.3434,0.4423,0.8465,0.1525,0.3697,0.132,0.3453,0.1296,0.0985,0.0827,0.1175,0.8051,1.3847,0.788,0.1961,0.1263,0.1633,0.1656,0.1254,0.0891,0.119,0.2616,0.4558,1.8698,0.3226,0.4071,0.337,0.2638,0.2377,0.162,0.1276,0.0851,0.081,0.7761,0.3557,0.3429,0.2652,0.2604,0.1528,0.1278,0.1181,0.086,0.1181,0.2912,0.2039,1.2034,0.6383,0.2852,0.3035,0.0835,0.2665,0.1342,0.0789,0.1353,1.2476,0.4096,0.5271,0.4088,0.1266,0.1281,0.1646,0.4593,0.8792,0.3281,0.1924,0.178,0.2001,0.1267,0.1832,0.1197,0.1501,0.2518,0.3123,0.0977,0.1079,0.4344,0.1601,0.2339,0.1474,0.1426]},
    midLower: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [0.1117,0.0856,0.0931,0.0596,0.0564,0.1336,0.3687,1.1051,0.278,0.1872,0.0962,0.071,0.1093,0.0857,0.0875,0.1761,0.1838,1.1563,0.3695,0.1513,0.093,0.094,0.1335,0.1203,0.0621,0.1851,0.3089,0.8839,0.0979,0.114,0.0889,0.1215,0.1648,0.0844,0.0839,0.1002,0.1965,0.4898,0.7613,0.1245,0.0883,0.0868,0.1158,0.1137,0.0927,0.2486,0.2425,0.1649,0.2395,0.1905,0.5145,0.159,0.0981,0.082,0.103,0.0994,0.0963,0.3615,0.3535,0.3693,0.3179,1.0634,1.9928,0.4673,0.0715,0.08,0.0613,0.1031,0.1132,0.0805,0.0899,0.6744,0.9174,0.1292,0.0816,0.1333,0.0687,0.1088,0.1043,0.1002,0.1107,0.4694,1.1355,2.0394,0.5423,0.1766,0.1171,0.0842,0.0789,0.0902,0.0836,0.1,0.1211,0.2559,0.58,1.6965,0.4325,0.1613,0.0789,0.1295,0.1325,0.1094,0.0699,0.122,0.0974,0.4983,3.1126,0.2198,0.1177,0.1335,0.1716,0.1131,0.0907,0.1603,0.2411,0.4709,0.2977,0.2005,0.1303,0.1049,0.1007,0.1043,0.1493,0.3545,0.3921,0.5761,0.0905,0.101,0.0874,0.126,0.0913,0.0857,0.0585,0.0662,0.7201,1.1889,0.4808,0.1255,0.098,0.1077,0.1757,0.111,0.0884,0.1358,0.2712,0.3901,1.3432,0.2222,0.1348,0.142,0.1059,0.1261,0.1279,0.111,0.092,0.0897,0.3965,0.1994,0.146,0.1084,0.1184,0.1301,0.1124,0.1225,0.0895,0.1028,0.2958,0.2078,0.5259,0.2573,0.1091,0.0835,0.0866,0.1441,0.1202,0.0684,0.1334,1.1705,0.3771,0.1917,0.143,0.1081,0.1215,0.1652,0.2641,0.5042,0.2086,0.0757,0.0906,0.0961,0.113,0.1093,0.0668,0.1535,0.2689,0.387,0.1046,0.0869,0.1826,0.1042,0.1031,0.1158,0.085]},
  },
  salinity: {
    surface: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [30.7205,31.1104,32.3467,32.6529,32.6732,33.1668,33.4012,33.2712,33.391,31.9405,29.303,29.7202,30.7651,31.8865,33.1227,33.2897,33.0144,33.1433,33.2797,32.4201,29.5885,29.9244,32.4755,32.6845,32.8628,33.1561,33.2099,33.2977,33.0433,31.9855,31.4589,30.0242,31.4696,31.8952,32.6384,32.5539,32.7581,33.0092,33.0494,33.0552,32.1076,31.1058,31.2266,32.3946,32.1907,32.8641,32.9406,33.2996,33.3969,33.4527,33.3761,32.8277,30.8036,28.0245,30.1967,30.9734,31.9851,31.6227,31.926,32.326,33.2371,33.1139,33.4486,33.6158,31.3821,28.0933,28.6978,30.9071,32.9427,33.1775,33.1516,32.9318,33.2834,32.4854,29.1075,28.518,29.6394,30.3584,31.8238,32.2368,32.5803,32.8112,33.0383,33.3285,33.2522,32.4129,31.924,30.9747,31.5545,31.9021,31.9873,32.6056,32.5352,33.004,33.0855,33.286,33.4332,30.0738,29.6548,30.4527,32.6076,32.4123,32.8152,33.0566,33.0979,33.4227,33.3189,32.837,31.6459,32.2838,32.2723,32.4406,32.8073,33.5253,33.4519,33.4493,33.1997,33.1033,30.3492,30.3249,31.2571,33.0391,33.1706,33.42,33.4831,33.335,33.3227,31.8828,29.6244,32.0227,31.8149,32.7164,32.5316,32.6135,33.2511,33.057,33.2764,31.4634,30.1159,30.467,32.7912,32.6856,33.1699,32.6964,33.1227,33.247,33.1808,33.1911,31.2535,29.8543,30.7267,31.0912,32.0064,32.0381,32.964,32.8917,33.4553,33.2854,31.806,30.8633,30.6685,30.8253,32.0616,32.5493,32.6382,32.7803,32.8749,33.0039,33.3671,33.2899,31.8342,28.6668,30.188,31.7355,32.7676,32.8704,32.8561,33.3824,32.9271,31.3657,31.1389,30.7767,31.4811,32.9946,33.0068,33.0828,33.11,30.6986,28.7332,28.0024,29.4006,31.9952,31.5887,32.5083,32.722,32.8251,32.6415,31.9273,31.0292,29.2482,30.3124,31.2854,31.3079]},
    mid: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [32.893,33.1852,33.1025,33.2206,33.6776,33.3616,33.4711,33.5388,33.4619,33.327,32.8835,32.612,32.8454,33.2334,33.2374,33.3034,33.611,33.4486,33.3979,33.2989,32.7388,32.963,32.8977,33.0862,33.5367,33.3656,33.3271,33.3289,33.3007,33.2073,32.7043,33.0469,33.0779,32.703,32.775,33.1353,33.3916,33.6103,33.5969,33.544,33.4279,33.234,33.0815,32.8977,33.2513,33.4352,33.8421,33.5051,33.4321,33.5004,33.4418,33.3763,33.1272,32.9172,32.6476,33.1039,32.7327,32.8158,32.9474,33.415,33.5609,33.7866,33.6915,33.6553,33.3202,33.0897,33.1591,33.3872,33.3297,33.1783,33.8071,33.6041,33.4969,33.3519,33.1536,32.6435,32.7034,32.7189,32.7212,32.9896,33.3737,33.4271,33.3905,33.3447,33.4536,33.3616,33.2103,32.7168,32.739,32.7495,32.7115,32.9638,32.9216,33.3678,33.5488,33.5309,33.4825,33.2459,32.9677,32.9503,32.9753,33.0186,33.4033,33.3606,33.6864,33.5125,33.4615,33.4235,33.2986,33.154,33.059,33.4301,33.7054,33.566,33.5459,33.5087,33.5003,33.3604,33.2087,33.2086,33.7768,33.4165,33.5339,33.4943,33.4888,33.5538,33.5142,33.2634,33.1944,33.1636,32.9435,32.6162,33.4944,33.5724,33.3128,33.4249,33.2999,33.3858,33.0569,33.1354,33.3297,33.3441,33.1745,33.1511,33.2727,33.356,33.3873,33.328,33.1587,32.9441,32.9214,33.2198,33.0987,33.1457,33.2361,33.5061,33.5952,33.4703,33.3883,33.1071,32.8358,33.0507,33.331,33.0178,32.9768,33.4346,33.4426,33.5697,33.5045,33.4747,33.0537,32.7692,33.3066,33.003,33.1928,32.9126,33.3818,33.4151,33.4428,33.2004,32.5651,33.1146,33.0537,33.0936,33.2178,33.1927,33.1825,33.0223,32.8866,31.5971,32.5452,33.0313,32.7406,33.0373,33.0593,33.2684,32.9465,32.6522,32.6932,31.4787,32.3628,32.5613,32.3199]},
    deep: {dates: ['2005-10','2005-11','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2008-02','2008-03','2008-05','2008-06','2008-07','2008-08','2008-10','2009-01','2009-04','2009-05','2009-06','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-07','2010-08','2010-09','2010-12','2011-01','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-12','2012-01','2012-04','2012-05','2012-06','2012-07','2012-09','2012-10','2012-12','2013-01','2013-03','2013-04','2013-05','2013-06','2013-09','2013-10','2013-11','2013-12','2014-01','2014-03','2014-04','2014-06','2014-07','2014-08','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-08','2019-09','2019-11','2020-01','2020-05','2020-06','2020-07','2020-08','2020-11','2021-01','2021-02','2021-03','2021-05','2021-06','2021-07','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-11','2023-12','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [33.4027,33.6806,33.7297,33.5413,33.5378,33.4282,33.2499,33.3391,33.6408,33.3499,33.3081,34.045,33.5434,33.4656,33.4511,33.327,34.1122,33.5408,33.3493,33.4091,33.364,33.2281,33.4594,32.9479,33.9338,33.9116,33.7203,33.444,33.313,33.5592,34.0175,34.2874,33.6295,33.4552,33.5526,33.5064,33.4227,33.4264,33.3134,33.594,33.4009,33.7445,34.0749,33.8177,33.721,33.7255,33.664,33.4308,33.2056,33.7374,33.6247,33.5403,33.4952,33.2933,33.274,33.256,33.3597,33.664,33.5635,33.359,33.5204,33.1465,33.3649,33.5449,32.9644,33.2039,33.5953,33.8466,33.5145,33.4873,33.4853,33.8914,33.5976,34.1463,33.5811,33.549,33.4632,33.5949,33.4671,33.3096,34.0028,33.5926,33.6249,33.5561,33.5843,33.4868,33.5249,33.3402,33.7078,33.5238,33.5166,33.6285,33.5609,33.4608,33.9706,33.7804,32.6054,34.1698,34.2409,33.4246,33.6311,33.3444,33.7774,33.4949,33.4313,33.7496,33.1911,33.2523,33.3865,33.5402,33.5731,33.3831,33.298,33.4505,33.5179,33.6056,33.6641,33.6329,33.6606,33.5189,33.7544,33.0475,34.1727,33.6778,33.6066,33.524,33.4104,32.7043,33.4269,32.9822,33.7947,33.4489,33.6127,33.3942,33.1665,33.7879,33.275,33.371,33.2493,33.2769,33.378,33.2562,33.5901,33.1763,33.715,33.1151,33.0397,33.3998,32.089,32.9509,33.1423,33.0066]},
    midUpper: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [31.7975,31.9754,32.4332,32.7004,32.8608,33.1801,33.4059,33.3628,33.413,32.8549,31.2059,31.0142,31.4757,32.2335,33.1471,33.2937,33.1781,33.2863,33.3337,33.0086,31.3565,31.9173,32.4949,32.7574,32.9527,33.1832,33.2326,33.3074,33.1786,32.7862,31.9575,31.9133,32.2729,32.2009,32.6574,32.6787,32.8222,33.1826,33.1177,33.2045,33.0858,32.4909,32.2689,32.4365,32.4634,32.95,33.1026,33.3395,33.4024,33.457,33.4017,33.0612,32.2325,31.1558,30.9362,32.0423,32.1959,32.1035,32.2195,32.4688,33.2594,33.2559,33.5433,33.6272,31.9723,31.588,31.5871,31.7297,33.0641,33.1692,33.2243,33.2487,33.3932,32.8469,31.9681,31.1174,31.0322,31.16,31.9847,32.5399,32.6839,32.9906,33.1789,33.3327,33.3139,33.1156,32.4373,31.8115,31.8663,32.039,32.26,32.6359,32.6669,33.046,33.1351,33.3258,33.4435,32.4939,31.473,32.0042,32.694,32.6381,32.8901,33.0741,33.1366,33.4414,33.3358,33.3343,32.7584,32.3988,32.602,32.632,32.9808,33.5354,33.4543,33.4624,33.4152,33.1807,32.1326,32.163,32.8886,33.1775,33.2892,33.4521,33.4831,33.4426,33.4312,32.6844,32.0745,32.1898,32.3456,32.6286,32.6499,32.7677,33.2666,33.1543,33.2807,32.849,32.3584,32.2378,32.9427,32.9636,33.1698,32.8843,33.151,33.27,33.2362,33.2484,32.5504,31.7313,31.9439,31.9583,32.4623,32.3165,32.986,32.985,33.5146,33.3251,32.7539,31.7945,32.1688,32.1123,32.4627,32.7213,32.8261,32.8523,33.0842,33.1098,33.4019,33.3935,32.1214,31.3701,32.1844,32.3937,32.8318,32.8722,32.9522,33.3848,33.1735,32.5703,31.3627,31.7436,32.403,33.0055,33.0473,33.0977,33.1131,32.1334,31.5733,29.9941,31.2414,32.1495,31.9912,32.7093,32.8196,32.8983,32.6815,31.9754,31.7589,30.4601,31.2947,31.6296,31.4244]},
    midCore: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [32.6797,32.9496,32.7676,33.0365,33.3077,33.2762,33.4282,33.447,33.4419,33.2785,32.6475,32.2031,32.613,33.0241,33.1983,33.3044,33.3405,33.4017,33.3747,33.2626,32.4663,32.9285,32.7501,32.976,33.2079,33.272,33.2994,33.3219,33.2577,33.1583,32.5217,32.8122,32.8847,32.5475,32.7338,32.9435,33.1164,33.3811,33.4073,33.4644,33.4113,33.2076,32.8834,32.6787,33.1453,33.181,33.532,33.4635,33.4212,33.482,33.4265,33.3605,32.992,32.7718,32.3143,33.0766,32.4543,32.4889,32.7432,33.0435,33.4625,33.6456,33.6445,33.6343,33.2388,32.8605,32.7827,32.92,33.2716,33.1623,33.4162,33.5315,33.4512,33.278,33.0548,32.8267,32.4627,32.6177,32.4727,32.7641,33.1626,33.2982,33.3137,33.3384,33.4002,33.3378,33.176,32.5549,32.4505,32.3011,32.5868,32.8409,32.8679,33.2115,33.3675,33.4456,33.4731,33.1995,32.8029,32.8889,32.8345,32.9078,33.1057,33.1837,33.376,33.4634,33.4125,33.4033,33.2082,33.0003,32.8734,33.1607,33.4636,33.5606,33.5038,33.4875,33.4821,33.2676,33.1384,33.0994,33.5415,33.3349,33.4153,33.4855,33.4845,33.5172,33.5055,33.2219,32.9929,32.8234,32.7376,32.6198,32.9754,33.0644,33.2849,33.321,33.2871,33.2198,33.0329,32.9832,33.2656,33.1846,33.1729,33.1115,33.2151,33.3097,33.3081,33.3095,33.109,32.8118,32.8378,33.0782,32.8946,32.8966,33.0406,33.257,33.5697,33.3964,33.2955,32.9583,32.885,32.7991,33.0635,32.9368,32.9533,33.0718,33.2975,33.4531,33.4515,33.4412,32.8797,32.72,33.1135,32.9116,33.1477,32.8828,33.1587,33.4007,33.3855,33.1289,32.1077,32.7971,33.0074,33.0386,33.1191,33.1702,33.144,32.9908,32.8062,31.8235,32.7251,32.7852,32.4916,32.8225,33.0538,33.0517,32.8873,32.4476,32.3513,31.2697,32.0707,32.2887,31.8114]},
    midLower: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [33.2519,33.5815,33.4694,33.4488,34.1506,33.4594,33.5133,33.6346,33.4849,33.4582,33.3731,33.2003,33.2715,33.5706,33.2813,33.3047,33.8761,33.5119,33.4487,33.3807,33.2,33.4289,33.0755,33.3805,33.8705,33.4638,33.3689,33.3378,33.3535,33.3234,33.197,33.4301,33.3659,32.935,32.8256,33.4115,33.7079,33.8471,33.8176,33.6643,33.5069,33.4109,33.3741,33.1344,33.4781,33.6999,34.1944,33.5655,33.4453,33.5213,33.46,33.4766,33.3945,33.3634,33.2083,33.3563,33.0295,33.2644,33.2272,33.927,33.686,33.9855,33.7522,33.6748,33.6419,33.5399,33.7439,34.0711,33.4209,33.1908,34.1894,33.7229,33.5493,33.5013,33.4545,33.2679,33.1949,33.0956,33.069,33.2283,33.7814,33.5992,33.4834,33.3513,33.5167,33.4472,33.3991,33.104,33.1741,33.1876,32.8839,33.1104,33.0821,33.5352,33.7511,33.6307,33.4966,33.426,33.3738,33.3729,33.1329,33.2015,33.7024,33.5346,34.0012,33.5592,33.5189,33.4546,33.4658,33.4059,33.2728,33.7671,34.0848,33.5756,33.592,33.5319,33.5292,33.4575,33.4691,33.5415,34.1485,33.52,33.6611,33.5084,33.4928,33.6002,33.5365,33.406,33.7188,33.5823,33.2167,32.6114,34.0059,34.0686,33.3405,33.5474,33.3121,33.6023,33.3417,33.4144,33.4744,33.5254,33.1765,33.2304,33.3351,33.4037,33.4698,33.356,33.3155,33.2727,33.3588,33.6883,33.3602,33.5605,33.4154,33.9157,33.628,33.5481,33.5758,33.4665,33.2776,33.4261,33.6806,33.1737,33.0223,33.8158,33.61,33.7717,33.5599,33.5131,33.3541,33.296,33.7672,33.1845,33.2943,32.9405,33.6149,33.4308,33.5342,33.3729,33.1069,33.5972,33.2567,33.1475,33.3171,33.2265,33.2218,33.2197,33.2007,33.0109,33.1013,33.3695,33.0544,33.2547,33.2001,33.4855,33.0383,32.922,33.1052,31.8195,32.7684,32.9269,32.8345]},
  },
  temperature: {
    surface: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [2.894,1.4172,1.6178,0.2247,0.2936,0.903,0.7623,1.1989,2.0916,4.8086,4.95,3.788,2.7661,1.2084,0.2917,0.0273,-0.2203,1.3688,1.9706,3.1957,4.1629,2.48,1.5336,0.1064,-0.197,-0.783,0.1092,0.3638,2.2891,3.8455,3.857,3.0655,1.5126,-0.298,0.3127,-0.2365,0.1318,1.0272,1.7128,2.6567,4.124,4.6009,3.6851,2.3233,1.0572,1.1943,0.1654,0.969,0.5173,0.9074,1.7023,3.3267,5.4185,5.147,4.6589,3.8484,2.4326,1.4245,0.5359,0.3889,1.1503,1.2852,1.7161,2.2919,4.903,5.3632,3.8283,2.3281,1.1227,-0.3235,-0.2919,0.8977,1.2975,3.3578,5.5885,3.7808,4.2469,2.7368,1.1832,0.4762,0.0384,0.846,1.3011,0.9197,2.7672,4.6294,4.2666,3.6735,3.2498,1.5873,0.4909,0.4565,-0.0817,0.4269,0.87,1.701,2.1658,5.329,5.7034,2.7395,1.8139,0.3357,-0.2582,-0.1402,-0.0992,-0.1361,0.5329,3.2247,3.7512,2.4888,1.9351,0.455,0.0611,-0.2533,0.0229,0.337,1.3521,2.4803,4.3445,2.8754,0.5984,1.3722,-0.0617,-0.0272,0.1751,1.2818,1.7731,4.2664,3.9196,2.3677,0.8882,0.6524,-0.9035,-0.0563,0.3914,1.2809,1.6513,3.6885,4.7454,3.3607,1.5124,0.4708,0.6367,-0.1985,0.3864,0.7056,1.5963,2.2177,4.2584,4.7866,3.9942,2.3016,1.6899,0.0571,0.61,0.5693,1.4382,2.3061,4.886,4.8592,3.0273,1.5777,1.209,0.8928,0.3557,0.0884,0.3843,0.9714,1.5717,2.6009,4.9467,5.9489,2.7713,2.2417,1.1281,0.0094,0.0388,0.8355,2.4778,4.0258,3.969,1.9417,1.0396,-0.0743,0.9796,1.4154,1.6689,3.9368,2.6389,3.5524,2.1881,1.6325,0.5121,-0.4841,0.0418,1.5316,2.1767,2.7406,2.9362,2.9597,2.6777,1.1962,0.0332]},
    mid: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [3.0334,3.4507,2.5562,1.3706,2.36,0.704,0.5703,0.772,1.9426,2.368,3.079,2.6574,2.8204,2.8175,0.2568,0.1404,0.896,0.802,1.6149,1.806,2.6972,3.4441,2.083,1.3146,1.491,-0.7965,-0.5017,0.2682,1.3039,1.6854,2.344,2.5165,2.4418,1.1239,0.6772,0.7903,0.7978,1.2161,1.2562,1.7493,1.8667,1.9868,2.4062,2.1997,2.7465,2.4243,2.6318,1.0696,0.533,0.9591,1.6199,2.099,2.6762,3.0419,3.9987,4.4588,3.069,3.1878,2.7005,2.8777,1.629,1.6696,1.3701,2.1533,2.8656,2.7644,3.108,3.7711,1.4606,-0.3066,1.7065,0.7021,1.0015,2.2274,2.3702,2.5431,2.5244,2.9127,2.2742,1.2742,1.5755,0.7437,0.9355,0.9127,2.1955,2.533,2.5411,2.6089,2.7354,2.2217,1.5061,0.9148,0.5151,0.6212,0.8623,1.3582,2.0375,2.3103,2.5669,3.1188,2.2477,1.4769,1.3655,0.0162,1.302,-0.1218,0.1742,1.158,1.8251,2.2875,2.1423,2.4584,2.4392,-0.3302,-0.2895,-0.1794,0.6982,1.7408,2.1088,2.7406,4.1055,2.2759,0.5927,-0.1576,0.1424,0.8717,1.0391,1.8205,2.9684,3.0557,2.187,0.7993,1.6865,1.912,0.4221,0.9237,1.59,2.6311,2.602,3.0613,2.6916,2.1017,0.6584,0.1232,0.1066,0.8153,1.2291,1.7362,1.9525,2.8693,3.4086,4.5223,3.4244,3.1106,1.1071,1.9036,1.1596,1.8838,2.1892,2.7382,2.879,3.6065,3.3482,1.7585,0.7799,1.7548,0.7752,1.1098,1.7159,2.0979,2.6042,2.7655,3.2712,2.7791,2.3721,0.2606,0.8696,0.8585,1.3476,1.9123,2.6704,3.4846,2.4877,-0.0172,0.5509,1.0881,1.4678,2.0048,2.0845,3.2966,2.7088,2.9882,1.8442,0.3628,0.4862,1.4879,1.6829,1.8272,2.6401,2.3577,2.0749,2.5296,1.5865]},
    deep: {dates: ['2005-10','2005-11','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2008-02','2008-03','2008-05','2008-06','2008-07','2008-08','2008-10','2009-01','2009-04','2009-05','2009-06','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-07','2010-08','2010-09','2010-12','2011-01','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-12','2012-01','2012-04','2012-05','2012-06','2012-07','2012-09','2012-10','2012-12','2013-01','2013-03','2013-04','2013-05','2013-06','2013-09','2013-10','2013-11','2013-12','2014-01','2014-03','2014-04','2014-06','2014-07','2014-08','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-08','2019-09','2019-11','2020-01','2020-05','2020-06','2020-07','2020-08','2020-11','2021-01','2021-02','2021-03','2021-05','2021-06','2021-07','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-11','2023-12','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [2.3443,3.9969,0.8526,1.6211,1.8,2.5985,2.3133,2.9261,3.2597,0.2669,0.0203,1.806,1.0383,1.2744,1.2741,2.1073,3.1898,-0.7872,0.2193,1.0801,1.2627,1.5683,2.7804,1.1918,1.5062,1.3466,1.655,2.1381,1.9609,3.2182,3.8513,3.5273,1.2421,0.5317,1.0328,1.7044,1.9703,2.8792,3.4235,4.4672,3.2869,1.8344,1.7625,1.1383,1.8995,2.3084,2.5787,1.5493,-0.3475,0.6538,1.2403,1.8267,1.8388,1.982,2.4867,3.1494,1.9477,0.6435,0.893,0.8951,2.1541,2.1884,2.7713,3.3193,1.9616,1.197,0.9026,1.0281,1.9401,1.5286,1.9643,2.5175,0.2566,2.3438,-0.0917,-0.0088,0.7992,1.6254,2.2674,2.1022,3.5545,-0.3782,-0.4857,-0.3963,0.7688,1.2393,1.8073,1.457,0.6022,-0.2128,0.1543,0.8471,0.7931,1.0212,3.4254,3.6776,1.2974,3.385,3.4573,0.4917,1.0999,1.4536,2.8891,3.0073,2.7101,3.1489,0.6891,0.2222,-0.2863,1.0132,1.1182,1.4059,2.0816,3.3788,4.3205,1.9103,1.0497,1.5766,1.9758,2.0196,4.0638,0.7993,3.3483,0.8676,1.6772,1.7569,1.8372,2.8504,2.8961,0.6724,1.6974,0.8748,1.1743,1.4035,1.8475,4.1246,0.051,0.4948,0.9649,1.32,1.9517,1.9667,3.8857,2.2562,2.2125,1.4037,1.3573,2.8278,2.4863,1.6626,3.1277,2.3892]},
    midUpper: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [3.3721,2.4151,1.7542,0.3027,0.7808,0.8189,0.7106,1.0401,2.1262,3.7138,4.2391,3.5114,2.7422,1.7305,0.276,0.0634,0.4723,0.9478,1.8275,2.5888,3.6752,3.0903,1.5509,0.3811,0.0304,-0.9023,-0.0063,0.2869,1.7351,2.631,3.4696,2.8615,1.9722,0.2077,0.3514,0.0893,0.2593,1.0553,1.4903,2.12,2.854,3.5427,3.2068,2.2988,1.5408,1.3989,0.9613,0.9405,0.5148,0.9084,1.6695,2.8903,4.2024,3.8719,4.6732,4.3776,2.7687,2.1473,1.432,1.0478,1.1908,1.6206,1.6181,2.2684,4.5307,3.8545,3.4606,2.7501,1.2393,-0.3123,0.194,0.6475,0.9776,3.0725,4.0844,3.0983,3.9151,3.0707,1.3773,0.8852,0.4416,0.9804,1.2535,0.9173,2.4498,3.4396,3.8441,3.4656,3.2264,1.5366,0.7992,0.4746,0.3642,0.4118,0.9882,1.5245,2.1043,3.9191,5.0242,3.1107,1.904,0.897,0.2714,-0.1763,0.1138,-0.1579,0.4851,1.6522,2.8739,2.4298,2.1177,0.9028,0.649,-0.2873,0.0225,0.215,0.8078,2.3583,3.4117,2.8547,2.7635,1.7246,0.4796,-0.1047,0.1543,0.9414,1.645,3.4767,3.2148,2.4599,1.5114,0.6684,-0.5352,0.2,0.4292,1.1306,1.637,2.8701,3.3863,3.389,1.98,1.2558,0.6455,-0.1242,0.4159,0.7246,1.3445,2.0494,3.2628,4.3393,4.1192,3.2584,2.4418,1.1767,0.6204,0.871,1.2779,2.1536,3.5346,4.5838,3.3584,2.7828,2.1323,1.2643,0.6793,0.4161,1.0709,0.8972,1.6541,2.357,4.4954,4.4872,2.9068,2.5586,1.3631,0.0355,0.1457,0.8333,1.87,3.058,3.9028,2.5588,1.9343,-0.0216,0.7249,1.3043,1.6525,3.7401,2.6854,3.6205,2.7007,1.8277,1.0879,0.0896,0.4485,1.0848,2.0502,2.6556,2.7392,2.7687,2.5753,1.6161,0.3219]},
    midCore: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [3.2424,3.1737,2.1965,1.1393,1.8038,0.6582,0.5638,0.7948,2.0465,2.5251,3.2462,2.8628,2.7712,2.613,0.2346,0.1784,0.5157,0.7764,1.7345,1.9877,2.9785,3.3673,1.8688,1.2232,0.5399,-0.7505,-0.3282,0.2795,1.3959,1.836,2.5999,2.5126,2.2442,0.7743,0.5518,0.438,0.7157,1.0306,1.2175,1.7407,1.9883,2.085,2.5524,2.2685,2.5908,1.8851,2.1488,1.0294,0.5287,0.9329,1.6111,2.13,3.0115,3.0374,4.356,4.5965,3.0618,2.6554,2.5794,2.4331,1.4967,1.775,1.4768,2.2457,3.0336,2.8385,2.9239,3.3113,1.4046,-0.2848,0.5587,0.7819,0.9549,2.4467,2.564,2.4757,2.7327,3.0496,1.9436,1.0467,1.2455,0.926,1.0076,0.9217,2.2648,2.7371,2.6633,2.7807,2.8239,1.642,1.2398,0.7687,0.5059,0.412,0.8901,1.3683,2.0697,2.5444,2.8764,3.1252,2.085,1.1737,0.6986,-0.1333,0.7154,-0.1503,0.2938,1.3383,1.8679,2.2858,2.1687,1.9524,1.8522,-0.3215,-0.1745,-0.0393,0.6968,2.1127,2.2429,2.6208,3.73,2.077,0.7436,-0.1451,0.1373,0.8762,1.2882,2.0719,2.731,2.7283,1.9398,0.6329,0.4748,0.7876,0.4151,0.8105,1.6241,2.4601,2.4398,3.0211,2.603,1.7312,0.6526,0.0845,0.3125,0.7711,1.2314,1.9202,2.2175,3.169,3.4967,4.2823,3.0496,2.6466,0.7038,1.4384,1.1964,2.0443,2.2258,2.999,2.9195,3.2774,2.8694,1.6445,0.8014,1.2389,0.8579,1.109,1.7112,2.214,2.9824,2.9423,2.974,2.7553,2.2893,0.0873,0.4835,0.8504,1.3301,2.2074,3.2396,3.1373,2.4464,-0.0129,0.5732,1.1079,1.5755,2.2056,2.1336,3.2544,2.6915,2.6146,1.5835,-0.0553,0.5742,1.175,1.7653,2.0775,2.6121,2.3779,2.2706,2.2681,0.8846]},
    midLower: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [2.8277,3.8398,2.9674,1.7408,3.1523,0.7119,0.5462,0.7036,1.8373,1.9819,2.7378,2.3517,2.8686,3.169,0.2678,0.1306,1.2321,0.7899,1.3996,1.5301,2.3115,3.67,2.3306,1.8379,2.4111,-0.806,-0.742,0.2569,1.1573,1.3978,1.6222,2.4502,2.666,1.597,0.825,1.263,0.9706,1.3709,1.2354,1.6814,1.59,1.5861,2.1503,2.1345,3.0892,2.9851,3.2837,1.1218,0.5395,0.9865,1.616,1.855,2.1509,2.88,3.6281,4.3761,3.1347,3.8886,3.0327,3.654,1.8037,1.6095,1.2502,2.0692,2.4236,2.4988,3.162,4.3144,1.5417,-0.32,2.7886,0.6601,1.0382,1.9143,1.9013,2.3177,2.1102,2.7907,2.7183,1.5022,2.2293,0.576,0.8245,0.9058,2.099,2.1154,2.1772,2.226,2.549,2.7417,1.8229,1.0991,0.5822,0.8013,0.8189,1.3174,2.0029,1.8357,1.8737,3.1138,2.4337,1.8648,2.0246,0.1534,1.9264,-0.0958,0.0333,0.9405,1.5884,2.2604,2.1296,3.1025,3.3681,-0.3444,-0.4276,-0.3506,0.6773,1.3718,1.7612,2.808,4.6818,2.5221,0.5152,-0.1764,0.1434,0.855,0.7538,1.3249,3.1279,3.3909,2.506,0.9354,2.9304,2.9968,0.4253,0.9576,1.5581,2.6969,2.5128,3.0228,2.9341,2.5151,0.6648,0.1979,-0.0912,0.8626,1.2047,1.5521,1.5079,2.3787,3.0624,5.0802,3.8678,3.9807,1.4708,2.6865,1.1118,1.7239,1.8977,2.1987,2.3485,4.0123,3.9068,1.9992,0.7857,2.4055,0.6617,1.161,1.7312,1.9695,1.978,1.9865,3.6247,2.8387,2.6273,0.4201,1.2692,0.8689,1.2555,1.4892,2.0486,3.8985,2.6619,-0.0192,0.5016,1.032,1.3597,1.527,1.9325,2.9913,2.7609,3.4662,2.1671,0.71,0.3466,1.7753,1.5554,1.4968,2.639,2.2627,1.8459,2.8842,2.3026]},
  },
  density: {
    surface: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [24.4786,24.8934,25.873,26.2008,26.2135,26.5773,26.7743,26.6427,26.6757,25.2704,23.1673,23.6081,24.5245,25.5298,26.5762,26.7242,26.5135,26.5277,26.5957,25.8073,23.4694,23.8715,25.9821,26.2321,26.3898,26.6512,26.6559,26.7127,26.3797,25.397,24.9837,23.9106,25.1759,25.6123,26.1848,26.1422,26.2903,26.443,26.43,26.3614,25.4731,24.6303,24.8117,25.8605,25.783,26.3159,26.435,26.68,26.785,26.807,26.692,26.1223,24.3037,22.139,23.9042,24.598,25.5247,25.3034,25.5991,25.9288,26.6185,26.5108,26.7497,26.84,24.8182,22.1726,22.7881,24.6702,26.3837,26.6499,26.6275,26.3885,26.6451,25.8467,22.9457,22.6556,23.5023,24.2021,25.4811,25.8527,26.1503,26.2946,26.449,26.7063,26.509,25.664,25.3138,24.6159,25.1135,25.5185,25.6508,26.1509,26.1204,26.4734,26.5139,26.6203,26.7037,23.7325,23.3661,24.2719,26.0686,26.0009,26.3543,26.5441,26.5754,26.8395,26.7213,26.1353,25.1416,25.7591,25.7913,26.018,26.3326,26.9278,26.8554,26.8371,26.5724,26.4148,24.0556,24.163,25.0503,26.4446,26.6319,26.8321,26.8729,26.6881,26.6447,25.2809,23.5181,25.5597,25.491,26.2295,26.1497,26.1819,26.6745,26.465,26.6163,25.0022,23.8316,24.2388,26.234,26.2137,26.5954,26.2556,26.5712,26.6536,26.5426,26.506,24.7822,23.6204,24.3884,24.8186,25.5935,25.7126,26.4312,26.3751,26.7746,26.5744,25.151,24.4123,24.4258,24.6524,25.6689,26.0811,26.1824,26.31,26.3711,26.4422,26.6946,26.5541,25.1718,22.5605,24.0559,25.3394,26.2426,26.3868,26.3732,26.7548,26.274,24.8931,24.7189,24.5893,25.214,26.4909,26.4441,26.4772,26.4816,24.372,22.9133,22.2623,23.4755,25.5902,25.3289,26.1154,26.2654,26.262,26.0693,25.4543,24.7218,23.2999,24.17,25.048,25.125]},
    mid: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [26.2011,26.3944,26.4066,26.589,26.8777,26.7461,26.842,26.8848,26.7442,26.6016,26.1891,26.0085,26.1815,26.4893,26.6708,26.7299,26.9336,26.8101,26.7165,26.6229,26.1054,26.2194,26.2812,26.4843,26.8272,26.8214,26.7779,26.7437,26.6597,26.5579,26.1046,26.3677,26.3977,26.1891,26.2753,26.5567,26.7642,26.9142,26.901,26.8243,26.7217,26.5564,26.4027,26.2732,26.5097,26.6805,26.9882,26.8393,26.8128,26.8426,26.7513,26.6632,26.4161,26.2198,25.9161,26.232,26.0701,26.1213,26.2707,26.6241,26.8462,27.0246,26.9693,26.8831,26.5548,26.3808,26.4069,26.5232,26.6727,26.6501,27.0263,26.9416,26.8364,26.6337,26.4621,26.043,26.0914,26.0727,26.1233,26.4114,26.6961,26.7965,26.7556,26.7202,26.7181,26.6163,26.4946,26.0954,26.1034,26.1494,26.1725,26.4133,26.4024,26.7556,26.8875,26.8412,26.7537,26.5406,26.2948,26.2391,26.3312,26.4206,26.7335,26.782,26.9653,26.9116,26.8556,26.7677,26.6208,26.4715,26.4067,26.6731,26.8915,26.9644,26.9461,26.9109,26.8582,26.6767,26.5278,26.4699,26.7997,26.6816,26.8913,26.8985,26.8794,26.8909,26.8479,26.5913,26.4467,26.4137,26.3096,26.1407,26.7734,26.8226,26.7228,26.7835,26.6399,26.6288,26.3679,26.3922,26.577,26.6351,26.5983,26.6079,26.7066,26.7352,26.7344,26.6518,26.4986,26.2542,26.1886,26.3138,26.3275,26.3881,26.6192,26.778,26.906,26.7553,26.6652,26.3945,26.168,26.2729,26.5188,26.4012,26.4322,26.7308,26.8072,26.8886,26.7951,26.7427,26.3619,26.1215,26.5084,26.3107,26.4946,26.4085,26.7491,26.7802,26.771,26.5354,25.9678,26.3347,26.3747,26.5684,26.6392,26.5871,26.5541,26.3843,26.273,25.1454,25.9508,26.3124,26.1722,26.5024,26.5149,26.6207,26.3499,26.1032,26.0749,25.1262,25.8543,25.9764,25.8481]},
    deep: {dates: ['2005-10','2005-11','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2008-02','2008-03','2008-05','2008-06','2008-07','2008-08','2008-10','2009-01','2009-04','2009-05','2009-06','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-07','2010-08','2010-09','2010-12','2011-01','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-12','2012-01','2012-04','2012-05','2012-06','2012-07','2012-09','2012-10','2012-12','2013-01','2013-03','2013-04','2013-05','2013-06','2013-09','2013-10','2013-11','2013-12','2014-01','2014-03','2014-04','2014-06','2014-07','2014-08','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-08','2019-09','2019-11','2020-01','2020-05','2020-06','2020-07','2020-08','2020-11','2021-01','2021-02','2021-03','2021-05','2021-06','2021-07','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-11','2023-12','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [26.6666,26.7391,27.034,26.832,26.8164,26.6664,26.5468,26.5674,26.7783,26.7613,26.74,27.2225,26.8725,26.7949,26.7833,26.6245,27.1609,26.9632,26.763,26.7621,26.7141,26.5845,26.676,26.3843,27.1551,27.1483,26.9732,26.7136,26.6243,26.7171,27.0203,27.2686,26.9284,26.8318,26.8804,26.7979,26.7115,26.6412,26.5023,26.6217,26.5845,26.9796,27.2497,27.0868,26.9559,26.9277,26.8566,26.7484,26.6741,27.0521,26.9249,26.8164,26.7794,26.6069,26.5523,26.4814,26.6627,26.9936,26.8979,26.7332,26.7752,26.4737,26.6013,26.6965,26.3449,26.5897,26.9228,27.117,26.7874,26.7951,26.7618,27.0434,26.9616,27.2617,26.9658,26.9358,26.8228,26.8734,26.7242,26.6109,27.0392,26.9883,27.0191,26.9596,26.922,26.8142,26.8054,26.6819,27.0313,26.9252,26.9016,26.9529,26.9016,26.8073,27.0259,26.8502,26.1025,27.1885,27.2383,26.8094,26.9382,26.6856,26.9206,26.6847,26.6595,26.8753,26.6103,26.6849,26.8176,26.8716,26.8914,26.7199,26.6032,26.6148,26.5766,26.8626,26.969,26.9086,26.9017,26.7849,26.7909,26.4883,27.1943,26.9914,26.8804,26.8085,26.7115,26.0672,26.6402,26.443,27.0259,26.8068,26.9197,26.7289,26.5153,26.8115,26.7118,26.766,26.6407,26.6403,26.6768,26.5784,26.6783,26.4923,26.927,26.5049,26.4475,26.6244,25.605,26.3556,26.3927,26.3462]},
    midUpper: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [25.2974,25.5183,25.9332,26.2353,26.3383,26.5931,26.7812,26.7266,26.6908,26.1075,24.7465,24.6615,25.0937,25.7742,26.5967,26.7258,26.6112,26.6704,26.6496,26.3298,24.9193,25.417,25.9965,26.2771,26.4522,26.6775,26.6798,26.7251,26.5316,26.1484,25.4158,25.4334,25.7893,25.8373,26.1981,26.2282,26.3357,26.5807,26.5002,26.5241,26.3694,25.834,25.6872,25.8961,25.9719,26.3718,26.5219,26.7139,26.7895,26.8105,26.715,26.3471,25.5653,24.7416,24.4892,25.3966,25.6665,25.6407,25.7828,26.0069,26.6339,26.6017,26.8327,26.8511,25.3255,25.0865,25.122,25.2958,26.4738,26.6428,26.663,26.6583,26.7538,26.1604,25.3667,24.7789,24.6392,24.8151,25.5984,26.0743,26.2144,26.431,26.5651,26.71,26.5859,26.3412,25.7632,25.3003,25.3647,25.6318,25.854,26.1744,26.2051,26.5081,26.5467,26.6649,26.7169,25.801,24.8776,25.4849,26.1315,26.1525,26.3897,26.5599,26.5965,26.8556,26.7376,26.6628,26.1066,25.8558,26.0422,26.1471,26.4425,26.9376,26.8574,26.8541,26.7829,26.4867,25.5608,25.6307,26.2192,26.5318,26.7003,26.8617,26.874,26.7969,26.741,25.9942,25.5319,25.6864,25.8792,26.158,26.232,26.2947,26.6849,26.5532,26.6208,26.1795,25.7425,25.6468,26.3245,26.3919,26.5949,26.404,26.5925,26.671,26.6052,26.5648,25.907,25.1536,25.3442,25.4351,25.9055,25.8757,26.4484,26.433,26.8332,26.6182,26.0444,25.179,25.5945,25.5986,25.9294,26.1969,26.3163,26.3516,26.5006,26.5318,26.7168,26.657,25.4473,24.8518,25.646,25.8415,26.279,26.387,26.446,26.7569,26.5181,25.9409,24.9032,25.322,25.896,26.4973,26.4917,26.4966,26.4853,25.5315,25.1762,23.8403,24.9098,25.7007,25.6212,26.2529,26.3237,26.3504,26.1108,25.4997,25.3201,24.281,24.9622,25.298,25.206]},
    midCore: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [26.0124,26.233,26.1683,26.4573,26.6305,26.68,26.8078,26.8094,26.7202,26.5511,25.986,25.6648,25.9997,26.3403,26.6403,26.7287,26.7398,26.7738,26.6895,26.5806,25.8646,26.1992,26.1791,26.4038,26.6306,26.7438,26.7484,26.7374,26.619,26.5083,25.9396,26.1801,26.2589,26.0862,26.2492,26.424,26.5478,26.742,26.7511,26.7609,26.6998,26.5289,26.2331,26.0924,26.439,26.5223,26.7832,26.8082,26.8041,26.8293,26.7393,26.6484,26.2807,26.1039,25.6154,26.1954,25.8484,25.9095,26.119,26.3709,26.7765,26.9036,26.9241,26.8588,26.4762,26.1913,26.1223,26.1967,26.6297,26.6361,26.797,26.8782,26.8022,26.5575,26.3685,26.1946,25.8827,25.98,25.9507,26.245,26.552,26.6818,26.6893,26.7144,26.6698,26.581,26.4577,25.9521,25.8654,25.8345,26.0904,26.3231,26.3597,26.6415,26.7397,26.7719,26.7435,26.4862,26.1396,26.1895,26.2309,26.3525,26.5399,26.6465,26.7562,26.8732,26.8099,26.7399,26.5457,26.3484,26.2558,26.5002,26.7481,26.9596,26.9068,26.8873,26.8435,26.5754,26.4616,26.3975,26.6526,26.6322,26.7869,26.8908,26.8761,26.861,26.8251,26.5408,26.3064,26.171,26.1639,26.153,26.4462,26.5006,26.7005,26.707,26.6271,26.51,26.3622,26.274,26.5339,26.5369,26.5971,26.5777,26.6498,26.7004,26.6704,26.6237,26.4398,26.1236,26.1141,26.2286,26.2003,26.2348,26.4877,26.6152,26.883,26.684,26.5887,26.2545,26.204,26.1037,26.3506,26.3447,26.4119,26.4797,26.6853,26.7949,26.7527,26.7067,26.193,26.0691,26.3816,26.2395,26.4657,26.3931,26.5945,26.7689,26.7262,26.4565,25.5558,26.1143,26.3411,26.5237,26.5583,26.5677,26.5157,26.3445,26.2046,25.3293,26.0958,26.1488,25.9914,26.3509,26.5056,26.4681,26.2965,25.9217,25.8037,24.9574,25.6059,25.7799,25.4871]},
    midLower: {dates: ['2005-10','2005-11','2005-12','2006-01','2006-02','2006-03','2006-04','2006-05','2006-06','2006-07','2006-08','2006-09','2006-10','2006-11','2007-02','2007-03','2007-04','2007-05','2007-06','2007-07','2007-09','2007-10','2007-11','2008-01','2008-02','2008-03','2008-04','2008-05','2008-06','2008-07','2008-08','2008-09','2008-10','2008-12','2009-01','2009-02','2009-03','2009-04','2009-05','2009-06','2009-07','2009-08','2009-09','2009-10','2009-11','2009-12','2010-01','2010-02','2010-03','2010-04','2010-05','2010-06','2010-07','2010-08','2010-09','2010-10','2010-11','2010-12','2011-01','2011-02','2011-03','2011-04','2011-05','2011-06','2011-07','2011-08','2011-09','2011-10','2011-12','2012-01','2012-03','2012-04','2012-05','2012-06','2012-07','2012-08','2012-09','2012-10','2012-12','2013-01','2013-02','2013-03','2013-04','2013-05','2013-06','2013-07','2013-08','2013-09','2013-10','2013-11','2013-12','2014-01','2014-02','2014-03','2014-04','2014-05','2014-06','2014-07','2014-08','2014-10','2014-11','2014-12','2015-01','2015-02','2015-03','2015-04','2015-05','2015-06','2015-08','2015-09','2015-10','2015-11','2015-12','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07','2016-09','2016-10','2016-12','2017-02','2017-03','2017-04','2017-05','2017-06','2017-08','2017-09','2017-10','2017-11','2018-01','2018-02','2018-03','2018-04','2018-05','2018-06','2018-07','2018-08','2018-09','2018-11','2018-12','2019-01','2019-02','2019-03','2019-04','2019-05','2019-06','2019-07','2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-05','2020-06','2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2022-02','2022-03','2022-05','2022-06','2022-07','2022-09','2022-11','2022-12','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-02','2024-03','2024-05','2024-06','2024-07','2024-08','2024-09','2024-10','2024-11','2024-12'], vals: [26.5056,26.6753,26.6671,26.7486,27.1923,26.8244,26.8775,26.9662,26.7708,26.7387,26.6102,26.5037,26.5179,26.73,26.7057,26.7316,27.126,26.8618,26.7726,26.7091,26.5059,26.5704,26.4054,26.6853,27.0319,26.9014,26.8218,26.7516,26.7122,26.6721,26.5554,26.6796,26.6104,26.3469,26.3079,26.751,27.0088,27.0944,27.08,26.926,26.8062,26.7295,26.6571,26.4678,26.6635,26.8465,27.2166,26.8847,26.8231,26.8579,26.7664,26.7628,26.6748,26.5903,26.3987,26.4417,26.3016,26.4168,26.468,26.9666,26.9345,27.1888,27.0264,26.9055,26.8511,26.7635,26.869,27.0159,26.7406,26.6608,27.2548,27.0398,26.8767,26.7783,26.7417,26.5602,26.5181,26.384,26.3678,26.5886,26.9778,26.9451,26.8373,26.7261,26.7763,26.7196,26.6765,26.4364,26.4666,26.4607,26.2901,26.5205,26.5282,26.8803,27.053,26.9244,26.7678,26.7236,26.6791,26.5769,26.4432,26.5413,26.93,26.9159,27.177,26.9481,26.9093,26.807,26.7727,26.6754,26.579,26.8922,27.1214,26.9729,26.9898,26.9379,26.8829,26.7815,26.7637,26.7251,27.0376,26.7454,26.9984,26.911,26.8827,26.9294,26.8842,26.7433,26.8513,26.7189,26.5044,26.1291,27.0977,27.1407,26.7451,26.8799,26.6522,26.7967,26.6037,26.6186,26.6725,26.7485,26.5997,26.6683,26.7669,26.7709,26.8025,26.6878,26.6582,26.5594,26.5689,26.6279,26.4956,26.6439,26.7401,27.047,26.9357,26.8298,26.8392,26.7287,26.5656,26.5344,26.7473,26.5094,26.4686,26.9896,26.9488,27.0479,26.8387,26.7835,26.6556,26.6083,26.8438,26.4511,26.5565,26.4229,26.9117,26.7923,26.851,26.7057,26.4523,26.6819,26.5235,26.6121,26.722,26.618,26.5932,26.5801,26.5361,26.299,26.3909,26.5424,26.4014,26.6596,26.6361,26.7754,26.4327,26.3433,26.4044,25.406,26.196,26.2414,26.2148]},
  },
};

// =================================================================
// OPTIONAL LIVE DATA OVERRIDE
// -----------------------------------------------------------------
// If you're running the separate monthly Python pipeline
// (refresh_nuuk_data.py + GitHub Actions), this pulls the live Earth
// Engine Table Asset it maintains and MERGES it into NUUK_DATA above -
// new months get appended, existing months get their values refreshed
// to the latest pull. If the asset doesn't exist (pipeline not set up
// yet, or this is the first run before it's ever executed), this fails
// silently and the tool just keeps using the hardcoded snapshot -
// nothing breaks either way. This is intentionally NOT a hard
// dependency.
// =================================================================
var LIVE_ASSET_ID = 'users/YOUR_USERNAME/nuuk_monthly_data'; // <-- must match EE_ASSET_ID in the pipeline's GitHub secret

function tryLoadLiveNuukData(onDone) {
  var fc = ee.FeatureCollection(LIVE_ASSET_ID);
  fc.evaluate(function(result, err) {
    if (err || !result || !result.features || result.features.length === 0) {
      onDone({loaded: false, reason: err ? String(err) : 'No live asset found (this is normal if the optional monthly pipeline is not set up).'});
      return;
    }
    var varKeys = Object.keys(VAR_LABELS);
    var depthKeys = Object.keys(DEPTH_LABELS);
    var touchedKeys = {};

    result.features.forEach(function(f) {
      var p = f.properties;
      var month = p.month;
      if (!month) return;
      varKeys.forEach(function(v) {
        depthKeys.forEach(function(d) {
          var key = v + '_' + d;
          if (p[key] === undefined || p[key] === null) return;
          var target = NUUK_DATA[v][d];
          var idx = target.dates.indexOf(month);
          if (idx === -1) {
            target.dates.push(month);
            target.vals.push(p[key]);
          } else {
            target.vals[idx] = p[key]; // refresh existing month with the latest pulled value
          }
          touchedKeys[key] = true;
        });
      });
    });

    // Re-sort every touched series chronologically - 'YYYY-MM' strings
    // sort correctly lexicographically, but new months were just
    // appended, not inserted in order.
    Object.keys(touchedKeys).forEach(function(key) {
      var parts = key.split('_');
      var d = parts.pop();
      var v = parts.join('_');
      var target = NUUK_DATA[v][d];
      var zipped = target.dates.map(function(date, i) { return [date, target.vals[i]]; });
      zipped.sort(function(a, b) { return a[0] < b[0] ? -1 : (a[0] > b[0] ? 1 : 0); });
      target.dates = zipped.map(function(z) { return z[0]; });
      target.vals = zipped.map(function(z) { return z[1]; });
    });

    onDone({loaded: true, monthsInAsset: result.features.length, seriesTouched: Object.keys(touchedKeys).length});
  });
}
// ---------------------------------------------------------------
// ---------------------------------------------------------------
// 1. STUDY SITE
// ---------------------------------------------------------------
print('[load 1/5] Study site setup starting...');
var NUUK = ee.Geometry.Point([-51.883, 64.117]);
var RING = NUUK.buffer(1000);
var REGION = NUUK.buffer(15000);
var ringFC = ee.FeatureCollection([ee.Feature(RING)]);
var siteFC = ee.FeatureCollection([ee.Feature(NUUK)]);

// Palette/vis definitions moved HERE (early) so the legend built below
// and the actual map layers built later both reference the SAME plain
// objects - they cannot drift out of sync, unlike the previous version
// where the legend hardcoded separate, independently-typed hex colors
// that had no structural connection to what visChl/visSst/visKd/
// visDepthClass actually specified for the map layers.
var visDepthClass = {min: 1, max: 3, palette: ['a6d9f0', '3182bd', '08306b']};
var visChl = {min: 0, max: 5, palette: ['08306b','2171b5','6baed6','c7e9b4','fee391','fe9929','d94801']};
var visSst = {min: -2, max: 12, palette: ['08306b','4292c6','fee391','fe9929','d94801']};
var visKd = {min: 0, max: 1, palette: ['f7fbff','9ecae1','3182bd','08306b']};
// min/max tightened from 0-3000 and gamma added specifically for this
// Arctic/glacial scene: snow and ice reflect very brightly (near-
// saturating a generic 0-3000 stretch to solid white) while deep water
// is very dark (crushed to solid black) - losing mid-tone detail in
// between. This stretch is tuned for Nuuk, not a universal default.
var visS2 = {bands: ['B4','B3','B2'], min: 0, max: 2200, gamma: 1.3};

// ---------------------------------------------------------------
// 2. UI HELPERS — ported directly from GeoMarineAnalysis's own
// MODULE D sidebar helpers (lbl/sHead/dynLbl/row/legRow/legDiv),
// same font sizes and conventions, so this tool matches its layout.
// ---------------------------------------------------------------
function lbl(txt, sz, col, bg, bld) {
  return ui.Label(txt, {fontSize: sz+'px', fontWeight: bld?'bold':'normal', color: col||'#111111', backgroundColor: bg||'rgba(0,0,0,0)', padding: '2px 4px', margin: '1px 0'});
}
function sHead(txt, bg, tc) { return lbl(txt, 10, tc||'#ffffff', bg||'#334455', true); }
function dynLbl(init, col) { return ui.Label(init, {fontSize: '9px', color: col||'#111111', backgroundColor: 'rgba(0,0,0,0)', padding: '1px 4px', margin: '0'}); }
function row(key, dl) { return ui.Panel([lbl(key+': ', 9, '#334466'), dl], ui.Panel.Layout.Flow('horizontal'), {margin: '1px 0'}); }
function legRow(hex, main, sub) {
  return ui.Panel([
    ui.Label('', {backgroundColor: hex, width: '16px', height: '12px', margin: '2px 5px 0 2px', padding: '0', border: '1px solid #888888'}),
    ui.Panel([lbl(main, 9, '#111111'), sub ? lbl(sub, 8, '#555555') : null].filter(Boolean), ui.Panel.Layout.Flow('vertical'), {margin: '0'})
  ], ui.Panel.Layout.Flow('horizontal'), {margin: '2px 0'});
}
function legDiv() { return ui.Label('', {margin: '3px 0 1px 0', backgroundColor: '#cccccc', height: '1px', stretch: 'horizontal'}); }
// Renders the REAL palette stops from a vis object ({min, max, palette})
// as a strip of colored swatches, with numeric min/max labels below - not
// a single arbitrary color chip. Takes the vis object directly so the
// legend is structurally guaranteed to match whatever the map layer
// actually uses (same object, not a separately-typed hex literal that
// can drift out of sync).
function gradientBar(title, vis, unit) {
  var n = vis.palette.length;
  var swatches = vis.palette.map(function(hex) {
    return ui.Label('', {backgroundColor: '#' + hex, height: '14px', stretch: 'horizontal', margin: '0', padding: '0'});
  });
  var bar = ui.Panel(swatches, ui.Panel.Layout.Flow('horizontal'), {margin: '2px 0 0 0', border: '1px solid #888888'});
  // A real number under EVERY swatch, not just the two ends - the exact
  // value that palette stop represents, linearly interpolated between
  // vis.min and vis.max the same way GEE itself maps values to colors.
  var valueLabels = vis.palette.map(function(hex, i) {
    var val = vis.min + (i / (n - 1)) * (vis.max - vis.min);
    var text = (Math.round(val * 10) / 10) + (unit || '');
    if (i === n - 1) text += '+';
    return lbl(text, 7, '#555555', null, false);
  });
  var ticks = ui.Panel(valueLabels, ui.Panel.Layout.Flow('horizontal'), {margin: '0'});
  return ui.Panel([lbl(title, 9, '#111111', null, true), bar, ticks], ui.Panel.Layout.Flow('vertical'), {margin: '3px 0'});
}
// Button style matching GeoMarineAnalysis's goToCoordsBtn convention:
// light tinted background + dark bold text + solid colored border.
// White text on a dark background rendered with poor contrast in
// testing. This is deliberately robust even in the worst case: if GEE
// fails to honor the 'color' override (a real, observed quirk), the
// browser's default dark text still reads fine against a LIGHT
// background - unlike the old white-on-dark case, where a failed
// override left near-invisible text.
function actionBtn(labelText, bg, tc, border, onClick) {
  return ui.Button({
    label: labelText,
    style: {fontSize: '12px', fontWeight: 'bold', margin: '3px 4px', backgroundColor: bg, color: tc,
      stretch: 'horizontal', padding: '8px 6px', border: '3px solid ' + border},
    onClick: onClick
  });
}

var STATUS = {SIG: 'significant', BORDER: 'borderline', NULL: 'not_significant'};
var COLORS = {}; COLORS[STATUS.SIG]='#2e7d32'; COLORS[STATUS.BORDER]='#f9a825'; COLORS[STATUS.NULL]='#757575';
var LABELS = {}; LABELS[STATUS.SIG]='SIGNIFICANT (p<0.05)'; LABELS[STATUS.BORDER]='BORDERLINE'; LABELS[STATUS.NULL]='NOT SIGNIFICANT';
var VAR_LABELS = {turbidity: 'Turbidity (FTU)', fluorescence: 'Fluorescence / Chlorophyll (ug/l)',
  salinity: 'Salinity', temperature: 'Temperature / SST proxy (C)', density: 'Density (kg/m3)'};
var DEPTH_LABELS = {surface: 'Surface (<20 dbar)', mid: 'Mid (20-300 dbar, coarse)', midUpper: 'Mid: 20-50 dbar', midCore: 'Mid: 50-150 dbar (best signal)', midLower: 'Mid: 150-300 dbar', deep: 'Deep (>300 dbar)'};

// Which variables have a REAL satellite counterpart, and which real
// dataset/band to fetch. Salinity and Density are deliberately absent -
// no satellite instrument measures either directly, so those stay
// in-situ-only, honestly, rather than faking a proxy.
var SATELLITE_MAP = {
  temperature: {datasetId: 'NASA/OCEANDATA/MODIS-Aqua/L3SMI', band: 'sst', tag: 'satellite: MODIS-Aqua SST'},
  fluorescence: {datasetId: 'NASA/OCEANDATA/MODIS-Aqua/L3SMI', band: 'chlor_a', tag: 'satellite: MODIS-Aqua chlorophyll-a'},
  turbidity: {datasetId: 'NASA/OCEANDATA/MODIS-Aqua/L3SMI', band: 'Kd_490', tag: 'satellite: MODIS-Aqua Kd_490 proxy'}
};

// Real, LIVE Earth Engine fetch of a monthly time series near the GF3
// station - one batched request (map + aggregate_array + a single
// evaluate() call), not 228 separate calls, following the same
// batching discipline your reference tool's own v10.107 fix used
// ("rebuilt from ~43 separate calls down to 3 batched calls").
function fetchSatelliteMonthlySeries(varKey, onDone) {
  var meta = SATELLITE_MAP[varKey];
  if (!meta) { onDone({error: 'No real satellite instrument measures ' + VAR_LABELS[varKey] + ' - switch this side to In-situ.'}); return; }
  print('[satellite] Starting live fetch for ' + varKey + ' (' + meta.datasetId + ', band=' + meta.band + '). This is a real Earth Engine round-trip, can take several seconds to a minute...');
  var tFetchStart = Date.now();
  var region = NUUK.buffer(2000);
  var startYear = 2005, endYear = 2024;
  var nMonths = (endYear - startYear + 1) * 12;
  var months = ee.List.sequence(0, nMonths - 1);
  var band = meta.band;
  var monthlyFC = ee.FeatureCollection(months.map(function(offset) {
    offset = ee.Number(offset);
    var year = ee.Number(startYear).add(offset.divide(12).floor());
    var month = offset.mod(12).add(1);
    var start = ee.Date.fromYMD(year, month, 1);
    var end = start.advance(1, 'month');
    // FIX (2nd attempt - first one didn't actually work): "Band pattern
    // 'Kd_490' did not match any bands" kept happening even after
    // wrapping .select() in ee.Algorithms.If(hasBand, ...). Turns out
    // that doesn't reliably protect against the error - Earth Engine
    // can validate/construct BOTH branches of an If() against the
    // image's band schema before picking one, so putting .select()
    // inside a conditional branch doesn't stop it from throwing if
    // that image lacks the band. Real fix: never call .select() on any
    // image that might lack the band at all. Flag each image with a
    // property using ONLY bandNames().contains() (which never throws,
    // regardless of whether the band exists), filter the collection
    // down to flagged=true first, and only THEN call .select() - by
    // that point every remaining image is guaranteed to have the band.
    var rawCol = ee.ImageCollection(meta.datasetId).filterDate(start, end).filterBounds(region);
    var flagged = rawCol.map(function(img) {
      return img.set('hasRequestedBand', img.bandNames().contains(band));
    });
    var withBand = flagged.filter(ee.Filter.eq('hasRequestedBand', true));
    var img = withBand.select([band]).mean();
    var meanVal = img.reduceRegion({reducer: ee.Reducer.mean(), geometry: region, scale: 4000, maxPixels: 1e9, bestEffort: true}).get(band);
    return ee.Feature(null, {date: start.format('YYYY-MM'), value: meanVal});
  }));
  var validFC = monthlyFC.filter(ee.Filter.notNull(['value']));
  var packaged = ee.Dictionary({dates: validFC.aggregate_array('date'), values: validFC.aggregate_array('value')});
  packaged.evaluate(function(result, err) {
    print('[satellite] Fetch for ' + varKey + ' returned after ' + (Date.now()-tFetchStart) + ' ms.');
    if (err) { onDone({error: 'Earth Engine fetch failed: ' + err}); return; }
    if (!result || !result.dates || result.dates.length < 12) {
      onDone({error: 'Live satellite fetch returned only ' + (result && result.dates ? result.dates.length : 0) + ' valid months - too sparse (cloud/ice coverage at this latitude, or this band genuinely missing from most scenes here). Try In-situ instead.'});
      return;
    }
    print('[satellite] ' + result.dates.length + ' valid real months fetched for ' + varKey + '.');
    onDone({dates: result.dates, vals: result.values});
  });
}

// Resolves either side's series: in-situ is instant (existing arrays),
// satellite is a real async Earth Engine call. onReady(series, error)
// is called either way, so the caller doesn't need to branch on sync
// vs async - it just waits for the callback.
function getSeriesForSide(varKey, depthKey, sourceKey, onReady) {
  if (sourceKey === 'satellite') {
    fetchSatelliteMonthlySeries(varKey, function(res) {
      if (res.error) { onReady(null, res.error); return; }
      onReady({dates: res.dates, vals: res.vals, tag: SATELLITE_MAP[varKey].tag}, null);
    });
  } else {
    var s = NUUK_DATA[varKey][depthKey];
    onReady({dates: s.dates, vals: s.vals, tag: 'in-situ CTD, ' + DEPTH_LABELS[depthKey]}, null);
  }
}

// ---------------------------------------------------------------
// 3. SIDEBAR — compact, matching GeoMarineAnalysis's 256-300px panel
// ---------------------------------------------------------------
print('[load 2/5] Sidebar construction starting...');
var panel = ui.Panel({style: {width: '300px', padding: '4px', backgroundColor: '#f0f5f0'}});

panel.add(lbl('Nuuk Deep-Water Diagnostic v5', 12, '#ffffff', '#1a4a2a', true));
panel.add(lbl('Godthaabsfjord (Nuup Kangerlua), SW Greenland - GEM MarineBasis GF3', 8, '#557755'));
panel.add(lbl('Everything below is computed LIVE in your browser on real nuuk_ctd.csv data - nothing is pre-computed.', 8, '#557755'));

// =================================================================
// 4. LIVE COUPLING ENGINE
// =================================================================
panel.add(sHead('COUPLING ENGINE - pick 2 real variables + depths', '#4a2a5a'));
panel.add(lbl('Preprocess -> lag-scan cross-correlation -> generate IAAFT surrogates -> evaluate p-value. Runs live, right here.', 8, '#5a3a6a'));

var varKeys = Object.keys(VAR_LABELS);
var depthKeys = Object.keys(DEPTH_LABELS);
var sourceItems = [{label: 'In-situ (real CTD sensor)', value: 'insitu'}, {label: 'Satellite (live fetch)', value: 'satellite'}];

panel.add(lbl('Variable A', 9, '#334466', null, true));
var selA = ui.Select({items: varKeys.map(function(k){return {label: VAR_LABELS[k], value: k};}), value: 'turbidity', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
var selADepth = ui.Select({items: depthKeys.map(function(k){return {label: DEPTH_LABELS[k], value: k};}), value: 'deep', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
var selASource = ui.Select({items: sourceItems, value: 'insitu', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(selA);
panel.add(selADepth);
panel.add(selASource);
panel.add(lbl('(Satellite ignores the depth dropdown above - it reads the sea surface. Only Temperature, Fluorescence, Turbidity have a real satellite equivalent; Salinity/Density do not.)', 8, '#888888'));

panel.add(lbl('Variable B', 9, '#334466', null, true));
var selB = ui.Select({items: varKeys.map(function(k){return {label: VAR_LABELS[k], value: k};}), value: 'fluorescence', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
var selBDepth = ui.Select({items: depthKeys.map(function(k){return {label: DEPTH_LABELS[k], value: k};}), value: 'deep', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
var selBSource = ui.Select({items: sourceItems, value: 'insitu', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(selB);
panel.add(selBDepth);
panel.add(selBSource);

panel.add(lbl('Surrogates (speed vs precision)', 9, '#334466', null, true));
var selSurr = ui.Select({items: [
  {label:'19 (fastest - min p=0.05)', value:19},
  {label:'20 (fast, default)', value:20},
  {label:'25', value:25},
  {label:'50', value:50},
  {label:'100', value:100},
  {label:'200', value:200},
  {label:'300 (slowest, highest precision)', value:300}
], value: 20, style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(selSurr);
panel.add(lbl('At N=19, the smallest possible p-value is 1/(19+1)=0.05 - the standard minimum surrogate count for a two-tailed test at alpha=0.05. Lower N = faster but coarser p-value steps.', 8, '#888888'));
panel.add(lbl('\u26A0 KEEP THIS BROWSER TAB IN THE FOREGROUND while the test runs. Chrome throttles background-tab timers (sometimes to once per MINUTE after a few minutes hidden) - switching tabs or minimizing mid-run is the single biggest cause of a run taking 10+ minutes instead of under 1.', 8, '#c62828'));

var runBtn = actionBtn('\u25B6 RUN COUPLING TEST', '#cce0ff', '#003388', '#0055cc', null);
panel.add(runBtn);

var stepStatus = dynLbl('Ready.', '#334466');
panel.add(stepStatus);

var engineResults = ui.Panel({style: {margin: '2px 0 0 0'}});
panel.add(engineResults);

function runFullPipeline(seriesA, tagA, seriesB, tagB, vA, dA, vB, dB, nSurr) {
  print('=== Coupling test started: ' + vA + '[' + tagA + '] x ' + vB + '[' + tagB + '], nSurr=' + nSurr + ' ===');
  stepStatus.setValue('Steps 2-4/4 running (lag-scan, IAAFT surrogates, p-value). This label will NOT update again until done - open the Console for live progress.');
  stepStatus.style().set('color', '#f9a825');

  // Yield again here, same reasoning as the outer button handler: this
  // status update needs a real paint before the heavy IAAFT loop starts.
  ui.util.setTimeout(function() {
    print('[coupling] Step 2 status painted. Running chunked fullCouplingTestChunked now (yields every ' +
      '10 surrogates)...');
    var tStart = Date.now();
    // FIX HISTORY: first found stepStatus.setValue() during every chunk
    // was triggering an expensive full-sidebar relayout (50 calls for
    // nSurr=100 -> ~87s just for 10% progress). Cut that to 3 total
    // widget updates for the whole run. BUT slowness persisted (a real
    // 50-surrogate run still took 10-15 minutes) - now suspecting
    // ui.util.setTimeout SCHEDULING ITSELF carries real per-call
    // overhead in this environment, independent of what runs inside the
    // callback. chunkSize raised from 2 to 10 (5x fewer setTimeout
    // calls for the same nSurr), and processChunk now logs BOTH the
    // actual compute time AND the gap between scheduling and firing for
    // every single chunk - check the Console for '[coupling] chunk
    // done: compute=Xms, setTimeout scheduling overhead=Yms' lines. If
    // Y is large while X stays small, that directly confirms setTimeout
    // itself is the bottleneck, not the math.
    try {
      // FIX (more aggressive): even ~4-5 widget updates during the loop
      // could still mean tens of seconds of relayout tax, given the
      // scale of the problem a real test exposed (87s for 10/100 with
      // 50 updates). Safest option: touch the widget ZERO times during
      // the loop. Progress lives entirely in the Console via print()
      // (confirmed cheap - plain text logging, not panel layout) until
      // the single final update in finishPipeline() below.
      fullCouplingTestChunked(
        seriesA.dates, seriesA.vals, seriesB.dates, seriesB.vals,
        {window: 24, nSurr: nSurr, chunkSize: 1},
        function onProgress(done, total) {
          print('[coupling] progress: ' + done + '/' + total + ' surrogates (' + (Date.now()-tStart) + ' ms elapsed)');
        },
        function onComplete(result) {
          var tElapsed = Date.now() - tStart;
          print('=== Coupling test finished in ' + tElapsed + ' ms ===');
          finishPipeline(result, tElapsed, vA, dA, vB, dB, tagA, tagB, nSurr);
        }
      );
    } catch (e) {
      stepStatus.setValue('ERROR: ' + e.message);
      stepStatus.style().set('color', '#c62828');
      print('=== Coupling test THREW an error: ' + e.message + ' ===');
    }
  }, 50);
}

function finishPipeline(result, tElapsed, vA, dA, vB, dB, tagA, tagB, nSurr) {
  if (result.error) {
    stepStatus.setValue('\u26A0 ' + result.error);
    stepStatus.style().set('color', '#c62828');
    return;
  }

  stepStatus.setValue('Done in ' + tElapsed + ' ms. Steps 1-4 complete (preprocess, lag-scan, ' + result.nSurr + ' IAAFT surrogates, p-value).');
  stepStatus.style().set('color', '#1a4a2a');

  var sig = result.pValue < 0.05;
  var status = sig ? STATUS.SIG : (result.pValue < 0.10 ? STATUS.BORDER : STATUS.NULL);

  engineResults.add(lbl(VAR_LABELS[vA] + ' [' + tagA + ']  x  ' + VAR_LABELS[vB] + ' [' + tagB + ']', 9, '#111111', null, true));
  engineResults.add(lbl(LABELS[status], 9, '#ffffff', COLORS[status], true));

  engineResults.add(row('Overlapping real months (n)', dynLbl(String(result.n))));
  engineResults.add(row('AC1 window', dynLbl(result.window + ' months')));
  engineResults.add(row('Best lag', dynLbl(result.bestLag + ' months')));
  engineResults.add(row('|corr| at best lag', dynLbl(result.maxAbs.toFixed(4), COLORS[status])));
  engineResults.add(row('IAAFT p-value', dynLbl(result.pValue.toFixed(4) + ' (' + result.nSurr + ' surrogates)', COLORS[status])));

  var lagLabels = result.lags.map(function(l){ return String(l); });
  var lagVals = result.lags.map(function(l){ return result.corrs[l]; });
  var lagChart = ui.Chart.array.values({array: ee.Array(lagVals.map(function(v){return [v];})), axis: 0, xLabels: ee.List(lagLabels)})
    .setChartType('ColumnChart')
    .setOptions({title: 'Cross-correlation vs lag (months)', legend: {position: 'none'}, colors: [COLORS[status]],
      hAxis: {title: 'Lag (months)'}, vAxis: {title: 'AC1-trajectory correlation'}, height: 160});
  engineResults.add(lagChart);

  if (result.surrStats && result.surrStats.length > 0) {
    var histChart = ui.Chart.array.values({array: ee.Array(result.surrStats.map(function(v){return [v];})), axis: 0})
      .setChartType('Histogram')
      .setOptions({title: 'IAAFT null distribution (' + result.nSurr + ' surrogates)', legend: {position: 'none'}, colors: ['#9575cd'],
        hAxis: {title: '|corr| under the null', viewWindow: {min: 0, max: 1}, ticks: [0, 0.2, 0.4, 0.6, 0.8, 1]},
        vAxis: {title: 'count'}, histogram: {bucketSize: 0.05, minValue: 0, maxValue: 1}, height: 160});
    engineResults.add(histChart);
    var pct = Math.round((1 - result.pValue) * 100);
    engineResults.add(lbl('Your real |corr|=' + result.maxAbs.toFixed(3) + ' exceeds ~' + pct + '% of the ' + result.nSurr +
      ' randomized (IAAFT) versions of this same data. ' + (sig ? 'Unusual enough to call significant.' : 'Not unusual enough to rule out chance.'), 8, '#555555'));
  }
}

runBtn.onClick(function() {
  var vA = selA.getValue(), dA = selADepth.getValue(), srcA = selASource.getValue();
  var vB = selB.getValue(), dB = selBDepth.getValue(), srcB = selBSource.getValue();
  var nSurr = selSurr.getValue();

  engineResults.widgets().reset([]);
  stepStatus.setValue('Step 1/4: fetching data (' + (srcA === 'satellite' || srcB === 'satellite' ? 'live satellite call in progress - can take 10-60s...' : 'in-situ, instant') + ')...');
  stepStatus.style().set('color', '#f9a825');

  // FIX: when both sides are in-situ, getSeriesForSide's callback fires
  // SYNCHRONOUSLY (no real async gap) - so the status update above and
  // the entire downstream computation (getSeriesForSide -> getSeriesForSide
  // -> runFullPipeline, which does all the IAAFT math AND builds every
  // result widget) previously ran as ONE uninterrupted call stack with
  // zero yields. Same root cause as the page-load hang: the "Step 1/4"
  // text set above never got a chance to paint before everything else
  // finished, matching the reported "no status" symptom exactly. This
  // setTimeout forces a real yield so the status text is visible first.
  ui.util.setTimeout(function() {
    print('[coupling] Status painted. Starting getSeriesForSide chain...');
    // Resolves side A first, then side B (each may be instant in-situ or
    // a real async satellite fetch) - only once BOTH are ready does the
    // actual (fast, synchronous) coupling math run.
    getSeriesForSide(vA, dA, srcA, function(seriesA, errA) {
      if (errA) {
        stepStatus.setValue('\u26A0 ' + errA);
        stepStatus.style().set('color', '#c62828');
        return;
      }
      getSeriesForSide(vB, dB, srcB, function(seriesB, errB) {
        if (errB) {
          stepStatus.setValue('\u26A0 ' + errB);
          stepStatus.style().set('color', '#c62828');
          return;
        }
        runFullPipeline(seriesA, seriesA.tag, seriesB, seriesB.tag, vA, dA, vB, dB, nSurr);
      });
    });
  }, 50);
});

// =================================================================
// 4b. MODEL FIT: AR(1) relaxation rate + BDS residual check
// -----------------------------------------------------------------
// Everything above (AC1, the Coupling Engine) is METRIC-based - a
// generic statistic computed from the raw data. This section is
// MODEL-based instead: fits an actual equation of motion to ONE real
// variable, x(t) = phi*x(t-1) + noise, and tracks the fitted
// parameter itself over time - giving a real relaxation rate in
// months, not an abstract 0-1 correlation number.
// =================================================================
panel.add(sHead('MODEL FIT: AR(1) relaxation rate + BDS check', '#5a3a1a'));
panel.add(lbl('IN PLAIN TERMS: this fits a simple physical model - "how strongly does this month get pulled back toward normal, based on how far last month was from normal?" A slower pull-back (declining relaxation rate) is the same critical-slowing-down idea as AC1, but expressed as a real number in months, not an abstract score. Then a second check (BDS) asks: after removing that simple pattern, is what is left over just random noise, or is there still real structure hiding in it? If the leftovers are genuinely random, the relaxation-rate number can be trusted. If not, the real dynamics are more complex than this simple model captures.', 8, '#5a3a1a'));

panel.add(lbl('Variable', 9, '#334466', null, true));
var mfVar = ui.Select({items: varKeys.map(function(k){return {label: VAR_LABELS[k], value: k};}), value: 'turbidity', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(mfVar);
panel.add(lbl('Depth', 9, '#334466', null, true));
var mfDepth = ui.Select({items: depthKeys.map(function(k){return {label: DEPTH_LABELS[k], value: k};}), value: 'deep', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(mfDepth);

panel.add(lbl('Fitting window (months per relaxation-rate estimate)', 9, '#334466', null, true));
var mfWindow = ui.Select({items: [
  {label: '24 (more windows, but each estimate is noisy)', value: 24},
  {label: '36 (validated default - reliable estimates, still enough windows)', value: 36},
  {label: '48 (most reliable estimates, fewer windows)', value: 48},
  {label: '60 (very reliable, few windows - only for long series)', value: 60}
], value: 36, style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(mfWindow);
panel.add(lbl('REAL TRADEOFF, tested directly: a shorter window gives more data points for the trend test, but each individual relaxation-rate estimate becomes much noisier (a 24-month estimate can easily be off by +/-0.18 around its true value). 36 months was the shortest window that reliably detected a genuine trend in direct testing - shorter windows risk missing a real signal entirely, not just losing some precision.', 8, '#888888'));

panel.add(lbl('Surrogates (speed vs precision)', 9, '#334466', null, true));
var mfSurr = ui.Select({items: [
  {label: '19 (fastest - min p=0.05)', value: 19},
  {label: '20 (fast, default)', value: 20},
  {label: '50', value: 50},
  {label: '100', value: 100}
], value: 20, style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(mfSurr);

var mfRunBtn = actionBtn('\u25B6 RUN MODEL FIT', '#ffe0b2', '#5a3a1a', '#ff9800', null);
panel.add(mfRunBtn);
var mfStatus = dynLbl('Ready.', '#5a3a1a');
panel.add(mfStatus);
var mfResults = ui.Panel({style: {margin: '2px 0 0 0'}});
panel.add(mfResults);

function finishModelFit(result, tElapsed, vKey, dKey, tag) {
  if (result.error) {
    mfStatus.setValue('\u26A0 ' + result.error);
    mfStatus.style().set('color', '#c62828');
    return;
  }
  mfStatus.setValue('Done in ' + tElapsed + ' ms.');
  mfStatus.style().set('color', '#1a4a2a');

  var tauSig = result.tauPValue < 0.05;
  var tauStatus = tauSig ? STATUS.SIG : (result.tauPValue < 0.10 ? STATUS.BORDER : STATUS.NULL);
  var bdsSig = result.bdsPValue < 0.05;
  var bdsStatus = bdsSig ? STATUS.SIG : (result.bdsPValue < 0.10 ? STATUS.BORDER : STATUS.NULL);

  mfResults.add(lbl(VAR_LABELS[vKey] + ' [' + tag + ']', 9, '#111111', null, true));
  mfResults.add(row('Overlapping real months (n)', dynLbl(String(result.n))));
  mfResults.add(row('Fitting window', dynLbl(result.window + ' months')));

  // --- PART 1: the relaxation-rate trend (the model-based CSD signal) ---
  mfResults.add(lbl('PART 1: Is the relaxation rate genuinely trending?', 9, '#5a3a1a', null, true));
  // FIX: a real, previously-manual caveat, now built into the tool
  // itself. IAAFT gives an honest p-value even at very few windows
  // (surrogates are drawn at the same n), but a naive read of
  // "SIGNIFICANT" alone can still mislead - a small handful of windows
  // makes a perfectly monotonic order (tau=+-1.0) far easier to hit by
  // real signal OR by small-sample luck than a longer trajectory would.
  // Directly demonstrated earlier: Turbidity Deep at window=24 gave
  // tau=+1.0000 from just 6 windows. Flag this explicitly, every time,
  // regardless of what the p-value says.
  if (result.realLambdas.filter(function(v) { return !isNaN(v); }).length < 10) {
    mfResults.add(lbl('\u26A0 CAUTION: only ' + result.realLambdas.filter(function(v) { return !isNaN(v); }).length + ' windows produced this trend estimate. With this few data points, a perfectly monotonic order (tau near +1 or -1) is much easier to hit by small-sample luck than with more windows - the p-value above already accounts for this honestly, but treat any "SIGNIFICANT" verdict here with extra skepticism. Consider a longer window or a longer real record before trusting this result on its own.', 8, '#c62828'));
  }
  mfResults.add(lbl(LABELS[tauStatus], 9, '#ffffff', COLORS[tauStatus], true));
  mfResults.add(row('AR(1) trend (Kendall tau)', dynLbl(result.realTau.toFixed(4), COLORS[tauStatus])));
  mfResults.add(row('Trend p-value', dynLbl(isNaN(result.tauPValue) ? 'n/a' : result.tauPValue.toFixed(4) + ' (' + result.tauNSurr + ' surrogates)', COLORS[tauStatus])));
  mfResults.add(lbl((result.realTau < 0 ? 'Negative tau means the relaxation rate is DECLINING over time - the system is taking longer to recover from small nudges. This is the real, model-based version of "critical slowing down." ' : 'Positive tau means the relaxation rate is INCREASING over time - the system is recovering FASTER from small nudges, the opposite of critical slowing down. ') +
    (tauSig ? 'This trend is unusual enough to call significant.' : 'This trend is not unusual enough to rule out chance at this surrogate count.'), 8, '#555555'));

  var lambdaLabels = result.realLambdas.map(function(v, i) { return 'W' + (i + 1); });
  var lambdaChart = ui.Chart.array.values({array: ee.Array(result.realLambdas.map(function(v){return [isNaN(v) ? 0 : v];})), axis: 0, xLabels: ee.List(lambdaLabels)})
    .setChartType('ColumnChart')
    .setOptions({title: 'Relaxation rate (lambda, per month) across non-overlapping windows', legend: {position: 'none'}, colors: [COLORS[tauStatus]],
      hAxis: {title: 'Window (chronological order)'}, vAxis: {title: 'lambda (per month) - lower = slower relaxation'}, height: 160});
  mfResults.add(lambdaChart);

  if (result.tauSurrStats && result.tauSurrStats.length > 0) {
    var tauHist = ui.Chart.array.values({array: ee.Array(result.tauSurrStats.map(function(v){return [v];})), axis: 0})
      .setChartType('Histogram')
      .setOptions({title: 'IAAFT null distribution of the trend statistic', legend: {position: 'none'}, colors: ['#9575cd'],
        hAxis: {title: 'Kendall tau under the null', viewWindow: {min: -1, max: 1}, ticks: [-1, -0.5, 0, 0.5, 1]},
        vAxis: {title: 'count'}, histogram: {bucketSize: 0.1, minValue: -1, maxValue: 1}, height: 160});
    mfResults.add(tauHist);
  }

  // --- PART 2: does the model actually capture what's going on? ---
  mfResults.add(lbl('PART 2: Did this simple model actually capture the real dynamics?', 9, '#5a3a1a', null, true));
  mfResults.add(lbl(bdsSig ? 'REAL STRUCTURE STILL LEFT OVER' : 'MODEL LOOKS ADEQUATE', 9, '#ffffff', COLORS[bdsStatus], true));
  mfResults.add(row('Fitted phi, OLS (whole series)', dynLbl(result.phi.toFixed(4))));
  mfResults.add(row('Fitted phi, robust (Theil-Sen)', dynLbl(result.phiRobust.toFixed(4))));
  if (Math.abs((result.phiRobust - result.phi) / result.phi) > 0.15) {
    mfResults.add(lbl('\u26A0 These differ by ' + (Math.abs((result.phiRobust - result.phi) / result.phi) * 100).toFixed(0) + '% - a few extreme months are meaningfully pulling the standard fit. The robust value is likely closer to "typical" month-to-month behavior.', 8, '#c62828'));
  }
  mfResults.add(row('BDS statistic on residuals', dynLbl(result.realBDS.toFixed(4), COLORS[bdsStatus])));
  mfResults.add(row('BDS p-value', dynLbl(isNaN(result.bdsPValue) ? 'n/a' : result.bdsPValue.toFixed(4) + ' (' + result.bdsNSurr + ' surrogates)', COLORS[bdsStatus])));
  mfResults.add(lbl(bdsSig ?
    'The leftover residuals still show structure this simple linear model did not capture. Treat the relaxation-rate number above as an incomplete picture, not the full story - the real dynamics here are more complex than a single-step relaxation model.' :
    'The leftover residuals look statistically like plain noise - this simple model is not obviously missing anything major. The relaxation-rate number above can be read at face value.', 8, '#555555'));

  if (result.bdsSurrStats && result.bdsSurrStats.length > 0) {
    var bdsMinV = result.bdsSurrStats[0], bdsMaxV = result.bdsSurrStats[0];
    for (var bi = 1; bi < result.bdsSurrStats.length; bi++) {
      if (result.bdsSurrStats[bi] < bdsMinV) bdsMinV = result.bdsSurrStats[bi];
      if (result.bdsSurrStats[bi] > bdsMaxV) bdsMaxV = result.bdsSurrStats[bi];
    }
    var bdsBinWidth = Math.max(0.001, (bdsMaxV - bdsMinV) / 15 || 0.01);
    var bdsHist = manualHistogramChart(result.bdsSurrStats, bdsBinWidth, 'IAAFT null distribution of the BDS statistic', '#9575cd', 'BDS statistic under the null');
    if (bdsHist) {
      mfResults.add(bdsHist);
    } else {
      mfResults.add(lbl('All ' + result.bdsSurrStats.length + ' surrogates produced the same BDS value - no variation to plot.', 8, '#555555'));
    }
  }

  // --- PART 3: did the relaxation rate genuinely shift between the
  // early and late halves of the record? A DIFFERENT question from
  // Part 1's continuous trend test - this asks specifically "is the
  // recent era different from the older era," reusing the SAME
  // already-computed lambda trajectory rather than re-fitting on a
  // thin subset. Uses full permutation, NOT IAAFT - this is a
  // question about temporal POSITION (early vs late), the same class
  // of question as exceedance clustering, and IAAFT's null hypothesis
  // would bias it the same way it did for clustering.
  var validLambdasForSplit = result.realLambdas.filter(function(v) { return !isNaN(v); });
  var elResult = permutationEarlyLateTest(validLambdasForSplit, 1000, Math.floor(Math.random() * 1e6));
  mfResults.add(lbl('PART 3: Has the relaxation rate shifted between the early and late halves of the record?', 9, '#5a3a1a', null, true));
  if (elResult.error) {
    mfResults.add(lbl(elResult.error, 8, '#888888'));
  } else {
    var elSig = elResult.pValue < 0.05;
    var elStatus = elSig ? STATUS.SIG : (elResult.pValue < 0.10 ? STATUS.BORDER : STATUS.NULL);
    mfResults.add(lbl(LABELS[elStatus], 9, '#ffffff', COLORS[elStatus], true));
    mfResults.add(row('Early half (n windows)', dynLbl(String(elResult.nEarly))));
    mfResults.add(row('Late half (n windows)', dynLbl(String(elResult.nLate))));
    mfResults.add(row('Late minus early (lambda)', dynLbl(elResult.realDiff.toFixed(4), COLORS[elStatus])));
    mfResults.add(row('Permutation p-value', dynLbl(elResult.pValue.toFixed(4) + ' (1000 reshuffles)', COLORS[elStatus])));
    mfResults.add(lbl((elResult.realDiff < 0 ?
      'The late half shows a LOWER relaxation rate than the early half - consistent with the system slowing down (getting closer to a tipping point) toward the end of the record, not just a gradual trend. ' :
      'The late half shows a HIGHER relaxation rate than the early half - the system looks like it is recovering FASTER recently, the opposite of critical slowing down. ') +
      (elSig ? 'This early-vs-late difference is unusual enough to call significant.' : 'This early-vs-late difference is not unusual enough to rule out chance - the early and late eras do not look reliably different at this window size.'), 8, '#555555'));
  }
}

mfRunBtn.onClick(function() {
  var vKey = mfVar.getValue(), dKey = mfDepth.getValue(), win = mfWindow.getValue(), nSurr = mfSurr.getValue();
  mfResults.widgets().reset([]);
  mfStatus.setValue('Step 1/2: fetching data...');
  mfStatus.style().set('color', '#f9a825');

  ui.util.setTimeout(function() {
    getSeriesForSide(vKey, dKey, 'insitu', function(series, err) {
      if (err) { mfStatus.setValue('\u26A0 ' + err); mfStatus.style().set('color', '#c62828'); return; }
      print('=== Model fit started: ' + vKey + '[' + series.tag + '], window=' + win + ', nSurr=' + nSurr + ' ===');
      mfStatus.setValue('Step 2/2: fitting AR(1), running BDS + IAAFT surrogates...');
      var tStart = Date.now();
      try {
        modelFitAnalysisChunked(series.dates, series.vals, {window: win, nSurr: nSurr, chunkSize: 1},
          function onProgress(done, total) {
            print('[model-fit] progress: ' + done + '/' + total + ' surrogates (' + (Date.now()-tStart) + ' ms elapsed)');
          },
          function onComplete(result) {
            var tElapsed = Date.now() - tStart;
            print('=== Model fit finished in ' + tElapsed + ' ms ===');
            finishModelFit(result, tElapsed, vKey, dKey, series.tag);
          }
        );
      } catch (e) {
        mfStatus.setValue('ERROR: ' + e.message);
        mfStatus.style().set('color', '#c62828');
        print('=== Model fit THREW an error: ' + e.message + ' ===');
      }
    });
  }, 50);
});

// =================================================================
// 4c. PEAKS-OVER-THRESHOLD: Generalized Pareto tail fit + exceedance
// clustering test
// -----------------------------------------------------------------
// A third, different kind of question from everything above: not "is
// there a trend" but "how extreme do the extreme events get, and do
// they cluster together in time." Fits a Generalized Pareto
// Distribution to values above a high threshold, and separately tests
// whether extreme months land next to each other more than chance
// reordering would produce.
// =================================================================
panel.add(sHead('PEAKS-OVER-THRESHOLD: extreme event analysis', '#7a2a2a'));
panel.add(lbl('IN PLAIN TERMS: everything above tests whether the middle of the data is behaving strangely. This asks a different question: how extreme do the WORST months actually get, and do bad months tend to cluster together (a real event lasting several months) or show up scattered and independent? A Generalized Pareto Distribution is the standard shape statisticians use to describe "how bad can it get" once you are already above a high threshold.', 8, '#7a2a2a'));
panel.add(lbl('HONEST LIMIT: Nuuk only has ~150-200 real months. At a 90th-percentile threshold that leaves roughly 15-20 extreme months to fit a distribution to - genuinely small. The tail-shape number (xi) below should be read as a rough indication, not a precise, confident figure; it can swing noticeably from run to run at this sample size, and that is a known property of this kind of estimate, not a flaw specific to this tool.', 8, '#888888'));

panel.add(lbl('Variable', 9, '#334466', null, true));
var evVar = ui.Select({items: varKeys.map(function(k){return {label: VAR_LABELS[k], value: k};}), value: 'turbidity', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(evVar);
panel.add(lbl('Depth', 9, '#334466', null, true));
var evDepth = ui.Select({items: depthKeys.map(function(k){return {label: DEPTH_LABELS[k], value: k};}), value: 'deep', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(evDepth);

panel.add(lbl('Threshold (percentile - lower = more exceedances, less extreme)', 9, '#334466', null, true));
var evPercentile = ui.Select({items: [
  {label: '85th percentile (more data, less extreme)', value: 85},
  {label: '90th percentile (default)', value: 90},
  {label: '95th percentile (fewer, more extreme events)', value: 95}
], value: 90, style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(evPercentile);

panel.add(lbl('Surrogates for the clustering test (full random reshuffles, not IAAFT - see Methodology)', 9, '#334466', null, true));
var evSurr = ui.Select({items: [
  {label: '19 (fastest - min p=0.05)', value: 19},
  {label: '20 (fast, default)', value: 20},
  {label: '50', value: 50},
  {label: '100', value: 100}
], value: 20, style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(evSurr);

var evRunBtn = actionBtn('\u25B6 RUN EXTREME VALUE ANALYSIS', '#ffcdd2', '#7a2a2a', '#e57373', null);
panel.add(evRunBtn);
var evStatus = dynLbl('Ready.', '#7a2a2a');
panel.add(evStatus);
var evResults = ui.Panel({style: {margin: '2px 0 0 0'}});
panel.add(evResults);

function finishExtremeValue(result, tElapsed, vKey, dKey, tag) {
  if (result.error) {
    evStatus.setValue('\u26A0 ' + result.error);
    evStatus.style().set('color', '#c62828');
    return;
  }
  evStatus.setValue('Done in ' + tElapsed + ' ms.');
  evStatus.style().set('color', '#1a4a2a');

  var clusterSig = result.clusterPValue < 0.05;
  var clusterStatus = clusterSig ? STATUS.SIG : (result.clusterPValue < 0.10 ? STATUS.BORDER : STATUS.NULL);

  evResults.add(lbl(VAR_LABELS[vKey] + ' [' + tag + ']', 9, '#111111', null, true));
  evResults.add(row('Overlapping real months (n)', dynLbl(String(result.n))));
  evResults.add(row('Threshold', dynLbl(result.percentile + 'th percentile (' + result.threshold.toFixed(3) + ')')));

  evResults.add(lbl('DESCRIPTIVE (not itself significance-tested - see Methodology for why)', 9, '#7a2a2a', null, true));
  evResults.add(row('Number of exceedances', dynLbl(String(result.nExceedances))));
  evResults.add(row('Mean excess over threshold', dynLbl(result.meanExcess.toFixed(4))));
  evResults.add(row('GPD shape (xi)', dynLbl(result.gpdXi.toFixed(3))));
  evResults.add(row('GPD scale (sigma)', dynLbl(result.gpdSigma.toFixed(3))));
  evResults.add(lbl((result.gpdXi > 0.05 ? 'Positive xi: a heavy tail - once you are past the threshold, how much further it goes has no natural ceiling. Extreme events here can be extreme in a real, open-ended way.' :
    result.gpdXi < -0.05 ? 'Negative xi: a bounded tail - there seems to be a natural ceiling on how far past the threshold this variable actually goes.' :
    'Xi near zero: an exponential-like tail - a fairly "ordinary" falloff, neither obviously bounded nor obviously open-ended.') + ' Remember the small-sample caveat above before reading too much into the exact number.', 8, '#555555'));

  evResults.add(lbl('DECLUSTERED comparison (adjacent exceedances merged into single events first)', 9, '#7a2a2a', null, true));
  evResults.add(lbl('Raw GPD theory assumes exceedances are roughly independent events. A sustained real event spanning 2-3 consecutive months otherwise gets counted as 2-3 separate "independent" observations, inflating the sample with duplicated information. This merges adjacent exceedances into one cluster first, keeping only each cluster\'s maximum, then refits - shown alongside the raw fit above, not in place of it, so you can see how much it actually moves.', 8, '#888888'));
  evResults.add(row('Real independent events (clusters)', dynLbl(String(result.nClusters) + ' (from ' + result.nExceedances + ' raw exceedances)')));
  if (isNaN(result.declGpdXi)) {
    evResults.add(lbl('Too few clusters (' + result.nClusters + ') for a reliable declustered GPD fit - need at least 5.', 8, '#888888'));
  } else {
    evResults.add(row('Declustered GPD shape (xi)', dynLbl(result.declGpdXi.toFixed(3))));
    evResults.add(row('Declustered GPD scale (sigma)', dynLbl(result.declGpdSigma.toFixed(3))));
    var xiShift = Math.abs(result.declGpdXi - result.gpdXi);
    evResults.add(lbl(xiShift > 0.15 ?
      'This moved xi by ' + xiShift.toFixed(3) + ' - a substantial shift. The raw fit above was likely distorted by real clustering (which the test below confirms or rules out) - the declustered numbers are the more defensible ones to actually read.' :
      'This moved xi by only ' + xiShift.toFixed(3) + ' - a modest shift. The raw and declustered pictures roughly agree here.', 8, '#555555'));
  }

  evResults.add(lbl('GENUINELY TESTED: do extreme months cluster together in time?', 9, '#7a2a2a', null, true));
  evResults.add(lbl(LABELS[clusterStatus], 9, '#ffffff', COLORS[clusterStatus], true));
  evResults.add(row('Adjacent-exceedance pairs', dynLbl(String(result.realPairCount) + ' of ' + result.totalExceed + ' exceedances')));
  evResults.add(row('Clustering p-value', dynLbl(isNaN(result.clusterPValue) ? 'n/a' : result.clusterPValue.toFixed(4) + ' (' + result.clusterNSurr + ' reshuffles)', COLORS[clusterStatus])));
  evResults.add(lbl(clusterSig ?
    'Extreme months land next to each other far more often than a random reshuffle of the same values would produce. This looks like real, sustained events - not scattered, independent bad months.' :
    'How often extreme months land next to each other is not unusual compared to a random reshuffle. No strong evidence of real clustering at this surrogate count - extreme months look more like scattered, independent events.', 8, '#555555'));

  if (result.clusterSurrStats && result.clusterSurrStats.length > 0) {
    var clusterHist = manualHistogramChart(result.clusterSurrStats, 1, 'Null distribution of adjacent-exceedance pairs (random reshuffles)', '#e57373', 'adjacent pairs under random reshuffling');
    if (clusterHist) {
      evResults.add(clusterHist);
    } else {
      var repeatedVal = result.clusterSurrStats[0];
      evResults.add(lbl('All ' + result.clusterSurrStats.length + ' random reshuffles produced exactly ' + repeatedVal + ' adjacent-exceedance pair(s) - no variation to plot.', 8, '#555555'));
    }
  }

  // Recency concentration: reuses the SAME exceedances/threshold as
  // clustering above, testing a different, complementary question -
  // are extreme months concentrated in the recent era, not just
  // whether they cluster together wherever they happen to fall.
  if (result.recency) {
    var rec = result.recency;
    var recSig = rec.pValue < 0.05;
    var recStatus = recSig ? STATUS.SIG : (rec.pValue < 0.10 ? STATUS.BORDER : STATUS.NULL);
    evResults.add(lbl('GENUINELY TESTED: are extreme months concentrated in the last 5 years?', 9, '#7a2a2a', null, true));
    evResults.add(lbl(LABELS[recStatus], 9, '#ffffff', COLORS[recStatus], true));
    evResults.add(row('Exceedances in last 5 years', dynLbl(String(rec.realRecentCount) + ' of ' + result.totalExceed + ' total, out of ' + rec.recentMonths + ' recent months')));
    evResults.add(row('Recency p-value', dynLbl(rec.pValue.toFixed(4) + ' (1000 reshuffles, one-sided)', COLORS[recStatus])));
    evResults.add(lbl(recSig ?
      'Extreme months are concentrated in the recent era far more than a random scatter across the whole record would produce. This is real evidence that things have gotten worse recently, not just an artifact of the recent period existing.' :
      'How many extreme months fall in the recent era is not unusual compared to a random scatter across the whole record. No strong evidence that the recent period is disproportionately extreme at this reshuffle count.', 8, '#555555'));
  }
}

evRunBtn.onClick(function() {
  var vKey = evVar.getValue(), dKey = evDepth.getValue(), pct = evPercentile.getValue(), nSurr = evSurr.getValue();
  evResults.widgets().reset([]);
  evStatus.setValue('Step 1/2: fetching data...');
  evStatus.style().set('color', '#f9a825');

  ui.util.setTimeout(function() {
    getSeriesForSide(vKey, dKey, 'insitu', function(series, err) {
      if (err) { evStatus.setValue('\u26A0 ' + err); evStatus.style().set('color', '#c62828'); return; }
      print('=== Extreme value analysis started: ' + vKey + '[' + series.tag + '], percentile=' + pct + ', nSurr=' + nSurr + ' ===');
      evStatus.setValue('Step 2/2: fitting GPD, testing clustering against random reshuffles...');
      var tStart = Date.now();
      try {
        extremeValueAnalysisChunked(series.dates, series.vals, {percentile: pct, nSurr: nSurr},
          function onProgress(done, total) {
            print('[extreme-value] progress: ' + done + '/' + total + ' reshuffles (' + (Date.now()-tStart) + ' ms elapsed)');
          },
          function onComplete(result) {
            var tElapsed = Date.now() - tStart;
            print('=== Extreme value analysis finished in ' + tElapsed + ' ms ===');
            finishExtremeValue(result, tElapsed, vKey, dKey, series.tag);
          }
        );
      } catch (e) {
        evStatus.setValue('ERROR: ' + e.message);
        evStatus.style().set('color', '#c62828');
        print('=== Extreme value analysis THREW an error: ' + e.message + ' ===');
      }
    });
  }, 50);
});

// =================================================================
// 4d. CROSS-VARIABLE EVENT COINCIDENCE ANALYSIS
// -----------------------------------------------------------------
// A genuinely different question from everything above: not "is
// there a trend in one variable" but "do extreme months in ONE
// variable coincide with extreme months in ANOTHER variable more
// than chance." Directly motivated by real evidence found in this
// data - Turbidity's extreme months overlap far more across depths
// than random scattering would predict (12 shared months between
// Surface and Mid vs ~1.8 expected by chance). This formalizes that
// kind of observation into a real, published methodology (Donges et
// al. 2016, Event Coincidence Analysis), adapted with permutation-
// based significance instead of the original paper's Poisson-
// asymptotic formula, for the same small-sample reasons used
// throughout this tool.
// =================================================================
function crossVariableECAAnalysis(seriesA, seriesB, opts) {
  opts = opts || {};
  var percentile = opts.percentile || 90;
  var tau = opts.tau !== undefined ? opts.tau : 1;
  var nPerm = opts.nPerm || 1000;
  var seed = opts.seed || Math.floor(Math.random() * 1e6);

  var aligned = alignSeries(seriesA.dates, seriesA.vals, seriesB.dates, seriesB.vals);
  if (aligned.dates.length < 30) {
    return {error: 'Only ' + aligned.dates.length + ' overlapping real months between these two series - need at least 30.'};
  }

  var detA = linearDetrend(deseasonalize(aligned.dates, aligned.a));
  var detB = linearDetrend(deseasonalize(aligned.dates, aligned.b));

  function toIndicator(det, pct) {
    var sorted = det.slice().sort(function(x, y) { return x - y; });
    var idx = Math.min(Math.floor((pct / 100) * sorted.length), sorted.length - 1);
    var thr = sorted[idx];
    return det.map(function(v) { return v > thr ? 1 : 0; });
  }
  var indA = toIndicator(detA, percentile);
  var indB = toIndicator(detB, percentile);

  // Both directions, since ECA is inherently directional (the
  // original paper's "precursor" vs "trigger" distinction) - A's
  // extremes coinciding with B is not necessarily the same rate as
  // B's extremes coinciding with A, since each has its own count of
  // extreme months.
  var aToB = eventCoincidenceTest(indA, indB, tau, nPerm, seed);
  var bToA = eventCoincidenceTest(indB, indA, tau, nPerm, seed + 999);

  return {n: aligned.dates.length, percentile: percentile, tau: tau, aToB: aToB, bToA: bToA};
}

// SIGNED version: same alignment/preprocessing as above, but splits
// Variable B's extremes into UP (top percentile) and DOWN (bottom
// percentile) SEPARATELY, then reuses the SAME already-validated
// eventCoincidenceTest against each - not new statistical machinery,
// just applying the existing test to two different definitions of "B's
// extreme." Distinguishes two real, competing, both-published
// mechanisms an unsigned coincidence test can't tell apart: does
// Variable A's extremes coincide with B spiking UP (consistent with a
// real upwelling/nutrient-pumping mechanism specifically documented at
// Godthaabsfjord) or B dropping DOWN (consistent with light-limitation
// from sediment loading)? Ported verbatim from the server-side version,
// same logic, tested there first.
function signedCrossVariableECAAnalysis(seriesA, seriesB, opts) {
  opts = opts || {};
  var percentile = opts.percentile || 90;
  var tau = opts.tau !== undefined ? opts.tau : 1;
  var nPerm = opts.nPerm || 1000;
  var seed = opts.seed || Math.floor(Math.random() * 1e6);

  var aligned = alignSeries(seriesA.dates, seriesA.vals, seriesB.dates, seriesB.vals);
  if (aligned.dates.length < 30) {
    return {error: 'Only ' + aligned.dates.length + ' overlapping real months between these two series - need at least 30.'};
  }

  var detA = linearDetrend(deseasonalize(aligned.dates, aligned.a));
  var detB = linearDetrend(deseasonalize(aligned.dates, aligned.b));

  function toIndicator(det, pct) {
    var sorted = det.slice().sort(function(x, y) { return x - y; });
    var idx = Math.min(Math.floor((pct / 100) * sorted.length), sorted.length - 1);
    var thr = sorted[idx];
    return det.map(function(v) { return v > thr ? 1 : 0; });
  }
  function toUpDownIndicators(det, pct) {
    var sorted = det.slice().sort(function(x, y) { return x - y; });
    var upIdx = Math.min(Math.floor((pct / 100) * sorted.length), sorted.length - 1);
    var downIdx = Math.max(0, Math.floor(((100 - pct) / 100) * sorted.length) - 1);
    var upThreshold = sorted[upIdx];
    var downThreshold = sorted[downIdx];
    return {
      up: det.map(function(v) { return v > upThreshold ? 1 : 0; }),
      down: det.map(function(v) { return v < downThreshold ? 1 : 0; }),
      upThreshold: upThreshold, downThreshold: downThreshold
    };
  }

  var indA = toIndicator(detA, percentile);
  var signedB = toUpDownIndicators(detB, percentile);

  var aToBUp = eventCoincidenceTest(indA, signedB.up, tau, nPerm, seed);
  var aToBDown = eventCoincidenceTest(indA, signedB.down, tau, nPerm, seed + 111);

  return {
    n: aligned.dates.length, percentile: percentile, tau: tau,
    upThreshold: signedB.upThreshold, downThreshold: signedB.downThreshold,
    aToBUp: aToBUp, aToBDown: aToBDown
  };
}

panel.add(sHead('CROSS-VARIABLE EVENT COINCIDENCE ANALYSIS', '#1a3a5a'));
panel.add(lbl('IN PLAIN TERMS: does an extreme month in one variable tend to show up around the same time as an extreme month in a DIFFERENT variable (or the same variable at a different depth)? This is a real, published method (Donges et al. 2016) adapted here with permutation-based testing, matched to how small this dataset actually is - the same approach used for clustering and recency elsewhere in this tool.', 8, '#1a3a5a'));

panel.add(lbl('Variable A (reference)', 9, '#334466', null, true));
var ecaVarA = ui.Select({items: varKeys.map(function(k){return {label: VAR_LABELS[k], value: k};}), value: 'turbidity', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(ecaVarA);
var ecaDepthA = ui.Select({items: depthKeys.map(function(k){return {label: DEPTH_LABELS[k], value: k};}), value: 'surface', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(ecaDepthA);

panel.add(lbl('Variable B (compared against)', 9, '#334466', null, true));
var ecaVarB = ui.Select({items: varKeys.map(function(k){return {label: VAR_LABELS[k], value: k};}), value: 'turbidity', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(ecaVarB);
var ecaDepthB = ui.Select({items: depthKeys.map(function(k){return {label: DEPTH_LABELS[k], value: k};}), value: 'mid', style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(ecaDepthB);

panel.add(lbl('Threshold (percentile)', 9, '#334466', null, true));
var ecaPercentile = ui.Select({items: [{label:'85th', value:85},{label:'90th (default)', value:90},{label:'95th', value:95}], value: 90, style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(ecaPercentile);

panel.add(lbl('Time tolerance (months)', 9, '#334466', null, true));
var ecaTau = ui.Select({items: [{label:'0 - exact same month only', value:0},{label:'1 - within +/-1 month (default)', value:1},{label:'2 - within +/-2 months', value:2}], value: 1, style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(ecaTau);

panel.add(lbl('Permutations', 9, '#334466', null, true));
var ecaPerm = ui.Select({items: [{label:'200 (fast)', value:200},{label:'1000 (default)', value:1000}], value: 1000, style: {stretch: 'horizontal', fontSize: '11px', margin: '1px 4px'}});
panel.add(ecaPerm);

var ecaRunBtn = actionBtn('\u25B6 RUN COINCIDENCE ANALYSIS', '#bbdefb', '#0d47a1', '#1976d2', null);
panel.add(ecaRunBtn);
var ecaStatus = dynLbl('Ready.', '#0d47a1');
panel.add(ecaStatus);
var ecaResults = ui.Panel({style: {margin: '2px 0 0 0'}});
panel.add(ecaResults);

function renderECADirection(title, res) {
  if (res.error) { ecaResults.add(lbl(title + ': ' + res.error, 8, '#888888')); return; }
  var sig = res.pValue < 0.05;
  var status = sig ? STATUS.SIG : (res.pValue < 0.10 ? STATUS.BORDER : STATUS.NULL);
  ecaResults.add(lbl(title, 9, '#0d47a1', null, true));
  ecaResults.add(lbl(LABELS[status], 9, '#ffffff', COLORS[status], true));
  ecaResults.add(row('Coincidence rate', dynLbl(res.nCoincide + ' of ' + res.nA + ' extreme months (' + (res.realRate*100).toFixed(0) + '%)', COLORS[status])));
  ecaResults.add(row('p-value', dynLbl(res.pValue.toFixed(4) + ' (' + res.permRates.length + ' permutations)', COLORS[status])));
}

function finishECA(result, tElapsed, vA, dA, tagA, vB, dB, tagB) {
  if (result.error) {
    ecaStatus.setValue('\u26A0 ' + result.error);
    ecaStatus.style().set('color', '#c62828');
    return;
  }
  ecaStatus.setValue('Done in ' + tElapsed + ' ms.');
  ecaStatus.style().set('color', '#1a4a2a');

  ecaResults.add(lbl(VAR_LABELS[vA] + ' [' + tagA + ']  vs  ' + VAR_LABELS[vB] + ' [' + tagB + ']', 9, '#111111', null, true));
  ecaResults.add(row('Overlapping real months (n)', dynLbl(String(result.n))));
  ecaResults.add(row('Threshold / tolerance', dynLbl(result.percentile + 'th percentile, +/-' + result.tau + ' month(s)')));

  renderECADirection('Does A\'s extremes coincide with B? (' + VAR_LABELS[vA] + ' -> ' + VAR_LABELS[vB] + ')', result.aToB);
  renderECADirection('Does B\'s extremes coincide with A? (' + VAR_LABELS[vB] + ' -> ' + VAR_LABELS[vA] + ')', result.bToA);

  var eitherSig = (result.aToB.pValue < 0.05) || (result.bToA.pValue < 0.05);
  ecaResults.add(lbl(eitherSig ?
    'At least one direction shows real coincidence - these two extreme-event patterns are not independent. Worth treating them as connected, not separate, coincidental phenomena.' :
    'Neither direction shows coincidence beyond what chance would produce at this sample size. These two extreme-event patterns look independent of each other here.', 8, '#555555'));
}

ecaRunBtn.onClick(function() {
  var vA = ecaVarA.getValue(), dA = ecaDepthA.getValue(), vB = ecaVarB.getValue(), dB = ecaDepthB.getValue();
  var pct = ecaPercentile.getValue(), tau = ecaTau.getValue(), nPerm = ecaPerm.getValue();
  ecaResults.widgets().reset([]);
  ecaStatus.setValue('Step 1/2: fetching both series...');
  ecaStatus.style().set('color', '#f9a825');

  ui.util.setTimeout(function() {
    getSeriesForSide(vA, dA, 'insitu', function(seriesA, errA) {
      if (errA) { ecaStatus.setValue('\u26A0 ' + errA); ecaStatus.style().set('color', '#c62828'); return; }
      getSeriesForSide(vB, dB, 'insitu', function(seriesB, errB) {
        if (errB) { ecaStatus.setValue('\u26A0 ' + errB); ecaStatus.style().set('color', '#c62828'); return; }
        print('=== ECA started: ' + vA + '[' + seriesA.tag + '] vs ' + vB + '[' + seriesB.tag + '], tau=' + tau + ', nPerm=' + nPerm + ' ===');
        ecaStatus.setValue('Step 2/2: testing coincidence against ' + nPerm + ' random reshuffles...');
        var tStart = Date.now();
        try {
          var result = crossVariableECAAnalysis(seriesA, seriesB, {percentile: pct, tau: tau, nPerm: nPerm});
          var tElapsed = Date.now() - tStart;
          print('=== ECA finished in ' + tElapsed + ' ms ===');
          finishECA(result, tElapsed, vA, dA, seriesA.tag, vB, dB, seriesB.tag);
        } catch (e) {
          ecaStatus.setValue('ERROR: ' + e.message);
          ecaStatus.style().set('color', '#c62828');
          print('=== ECA THREW an error: ' + e.message + ' ===');
        }
      });
    });
  }, 50);
});

// =================================================================
// 4e. SIGNED COINCIDENCE: UP vs DOWN
// -----------------------------------------------------------------
// Uses the same A/B/Depth/Permutations selections above. Splits
// Variable B's extremes into UP (spike) and DOWN (crash) separately -
// distinguishes an upwelling/nutrient-pumping mechanism (A coincides
// with B spiking up) from a light-limitation mechanism (A coincides
// with B dropping down), which the unsigned coincidence test above
// cannot tell apart.
// =================================================================
panel.add(sHead('SIGNED COINCIDENCE: UP vs DOWN', '#0d5c46'));
panel.add(lbl('Uses the same A/B/Depth/Permutations selections above. Splits Variable B\'s extremes into UP (spike) and DOWN (crash) separately - an upwelling/nutrient-pumping mechanism would show up as A coinciding with B spiking UP; a light-limitation mechanism would show up as A coinciding with B dropping DOWN. The unsigned test above cannot tell these apart.', 8, '#0d5c46'));
var signedRunBtn = actionBtn('\u25B6 RUN SIGNED COINCIDENCE', '#b2dfdb', '#0d5c46', '#00897b', null);
panel.add(signedRunBtn);
var signedStatus = dynLbl('Ready.', '#0d5c46');
panel.add(signedStatus);
var signedResults = ui.Panel({style: {margin: '2px 0 0 0'}});
panel.add(signedResults);

function renderSignedDirection(title, res) {
  if (res.error) { signedResults.add(lbl(title + ': ' + res.error, 8, '#888888')); return; }
  var sig = res.pValue < 0.05;
  var status = sig ? STATUS.SIG : (res.pValue < 0.10 ? STATUS.BORDER : STATUS.NULL);
  signedResults.add(lbl(title, 9, '#0d5c46', null, true));
  signedResults.add(lbl(LABELS[status], 9, '#ffffff', COLORS[status], true));
  signedResults.add(row('Coincidence', dynLbl(res.nCoincide + ' of ' + res.nA + ' extreme months (' + (res.realRate*100).toFixed(0) + '%)', COLORS[status])));
  signedResults.add(row('p-value', dynLbl(res.pValue.toFixed(4) + ' (' + res.permRates.length + ' permutations)', COLORS[status])));
}

function finishSignedECA(result, tElapsed, vA, dA, tagA, vB, dB, tagB) {
  if (result.error) {
    signedStatus.setValue('\u26A0 ' + result.error);
    signedStatus.style().set('color', '#c62828');
    return;
  }
  signedStatus.setValue('Done in ' + tElapsed + ' ms.');
  signedStatus.style().set('color', '#1a4a2a');

  signedResults.add(lbl(VAR_LABELS[vA] + ' [' + tagA + ']  vs  ' + VAR_LABELS[vB] + ' [' + tagB + ']', 9, '#111111', null, true));
  signedResults.add(row('UP threshold', dynLbl(result.upThreshold.toFixed(3) + ' (deseasonalized/detrended anomaly units)')));
  signedResults.add(row('DOWN threshold', dynLbl(result.downThreshold.toFixed(3))));

  renderSignedDirection('A coincides with B spiking UP (upwelling-consistent)', result.aToBUp);
  renderSignedDirection('A coincides with B dropping DOWN (light-limitation-consistent)', result.aToBDown);

  var upSig = result.aToBUp.pValue < 0.10;
  var downSig = result.aToBDown.pValue < 0.10;
  signedResults.add(lbl(
    (upSig && !downSig) ? 'The signal is directional: coincidence with B spiking UP, not with B dropping DOWN. Consistent with an upwelling/nutrient-pumping mechanism, not light-limitation.' :
    (downSig && !upSig) ? 'The signal is directional: coincidence with B dropping DOWN, not with B spiking UP. Consistent with light-limitation, not upwelling.' :
    (upSig && downSig) ? 'Both directions show some coincidence - possibly both mechanisms operating, or the split is not cleanly separating them at this sample size.' :
    'Neither direction shows coincidence beyond what chance would produce - no support for either mechanism from this specific test.', 8, '#555555'));
}

signedRunBtn.onClick(function() {
  var vA = ecaVarA.getValue(), dA = ecaDepthA.getValue(), vB = ecaVarB.getValue(), dB = ecaDepthB.getValue();
  var pct = ecaPercentile.getValue(), tau = ecaTau.getValue(), nPerm = ecaPerm.getValue();
  signedResults.widgets().reset([]);
  signedStatus.setValue('Step 1/2: fetching both series...');
  signedStatus.style().set('color', '#f9a825');

  ui.util.setTimeout(function() {
    getSeriesForSide(vA, dA, 'insitu', function(seriesA, errA) {
      if (errA) { signedStatus.setValue('\u26A0 ' + errA); signedStatus.style().set('color', '#c62828'); return; }
      getSeriesForSide(vB, dB, 'insitu', function(seriesB, errB) {
        if (errB) { signedStatus.setValue('\u26A0 ' + errB); signedStatus.style().set('color', '#c62828'); return; }
        print('=== Signed ECA started: ' + vA + '[' + seriesA.tag + '] vs ' + vB + '[' + seriesB.tag + '], tau=' + tau + ', nPerm=' + nPerm + ' ===');
        signedStatus.setValue('Step 2/2: testing UP/DOWN coincidence against ' + nPerm + ' random reshuffles...');
        var tStart = Date.now();
        try {
          var result = signedCrossVariableECAAnalysis(seriesA, seriesB, {percentile: pct, tau: tau, nPerm: nPerm});
          var tElapsed = Date.now() - tStart;
          print('=== Signed ECA finished in ' + tElapsed + ' ms ===');
          finishSignedECA(result, tElapsed, vA, dA, seriesA.tag, vB, dB, seriesB.tag);
        } catch (e) {
          signedStatus.setValue('ERROR: ' + e.message);
          signedStatus.style().set('color', '#c62828');
          print('=== Signed ECA THREW an error: ' + e.message + ' ===');
        }
      });
    });
  }, 50);
});

// =================================================================
// 5. METHODOLOGY — Definition / Data Source & Collection & Real
// Coverage / CSD-IAAFT Methodology Applied.
// FIX: the previous "collapsed panel, shown:false" approach did NOT
// avoid the layout cost - it only DEFERRED it to whenever the toggle
// was clicked (confirmed: that toggle click itself is what hung the
// page next). GEE's ui.Panel appears to do expensive layout work when
// ANY widget's `shown` state flips, regardless of how long it was
// hidden first. The actual fix: don't build ANY hidden widget tree at
// all. This button prints the methodology straight to the Console
// instead - print() is plain text logging, not panel layout, so it
// carries none of that risk.
// =================================================================
var methodBtn = actionBtn('\u2139 METHODOLOGY (click - prints to Console below)', '#b2dfdb', '#00695c', '#00897b', function() {
  print('--- METHODOLOGY: 1. Definition ---');
  print('Turbidity = water cloudiness from suspended particles (FTU), a proxy for glacial sediment. Fluorescence = chlorophyll/phytoplankton proxy. AC1 (lag-1 autocorrelation) = how strongly this month resembles last month; rising AC1 is the "critical slowing down" signature.');
  print('--- 2. Data Source, Collection, Real Coverage ---');
  print('GEM MarineBasis Nuuk, real in-situ CTD (SeaBird 19Plus), station GF3. api.g-e-m.dk, DOI:10.17897/KMEK-TK21. Real casts 2005-10-05 to 2024-12-11 (~19.2 yrs), depth 1-392 dbar. Monthly coverage: Surface/Mid ~195-203 months, Deep (>300 dbar) ~147-152 months.');
  print('--- 3. CSD/IAAFT Methodology ---');
  print('Preprocess (align, deseasonalize from first-half climatology, detrend) -> Lag-scan (rolling 24mo AC1 trajectory per series, correlated at lags -6..+6, keep max |corr|) -> Generate IAAFT surrogates (phase-randomized fakes sharing the real spectrum) -> Compare & Evaluate (p = fraction of surrogates meeting/beating the real stat).');
});
panel.add(methodBtn);


// =================================================================
// 6. LEGEND — merged from the old floating map overlay per feedback
// (it was blocking the map view); now lives in the scrollable sidebar.
// =================================================================
panel.add(sHead('LEGEND', '#37474f'));
panel.add(legRow(COLORS[STATUS.SIG], LABELS[STATUS.SIG], 'p < 0.05'));
panel.add(legRow(COLORS[STATUS.BORDER], LABELS[STATUS.BORDER], '0.05 <= p < 0.10'));
panel.add(legRow(COLORS[STATUS.NULL], LABELS[STATUS.NULL], 'p >= 0.10'));
panel.add(legDiv());
panel.add(legRow('#FFFFFF', 'Study zone border [in-situ]', '1km ring, white outline drawn on the map, centered on the real GF3 CTD station'));
panel.add(legRow('#2e7d32', 'GF3 station marker [in-situ]', 'real CTD cast location'));
panel.add(legDiv());
panel.add(lbl('Depth layer colors [satellite: NOAA ETOPO1] (matches Surface/Mid/Deep everywhere in this tool):', 8, '#333333'));
panel.add(legRow('#' + visDepthClass.palette[0], 'Surface', '< 20 dbar'));
panel.add(legRow('#' + visDepthClass.palette[1], 'Mid', '20-300 dbar'));
panel.add(legRow('#' + visDepthClass.palette[2], 'Deep', '> 300 dbar'));
panel.add(legDiv());
panel.add(lbl('Satellite layer color scales (enable a layer in native Layers panel to see it). These bars use the EXACT same palette as the map layer - not a separate approximation:', 8, '#333333'));
panel.add(gradientBar('Chlorophyll-a [MODIS-Aqua]', visChl, ' mg/m3'));
panel.add(gradientBar('SST [MODIS-Aqua]', visSst, 'C'));
panel.add(gradientBar('Turbidity proxy Kd_490 [MODIS-Aqua]', visKd, ''));
panel.add(lbl('KNOWN LIMIT: ETOPO1 is ~1.85km/pixel - can under-resolve this narrow fjord\'s true depth right at the station point. Trust the real CTD depth (392 dbar, in the map info box) over this coarse public raster.', 8, '#c62828'));
panel.add(lbl('p > 0.05 means "not distinguishable from noise" - never "confirmed absent."', 8, '#888888'));

// =================================================================
// 6b. WHAT THE SATELLITE LAYERS ACTUALLY ARE — same fix as METHODOLOGY
// above: prints to Console instead of building any hidden panel, since
// toggling `shown` on a hidden panel is what caused the hang, not just
// having long text visible.
// =================================================================
var satBtn = actionBtn('\uD83D\uDEF0 WHAT THE SATELLITE LAYERS SHOW (click - prints to Console below)', '#d7ccc8', '#4e342e', '#8d6e63', function() {
  print('--- Sentinel-2 ---');
  print('A real optical photo, ~10m resolution - shows actual coastline, ice, water color, like a zoomed-out aerial photo.');
  print('--- Black or transparent patches on Sentinel-2 ---');
  print('If any small dark/no-data patches remain, that means every single Sentinel-2 pass over that exact spot in the whole date window was cloudy or missing - a real gap in what the satellite could see, not a rendering bug. Widened the cloud filter (<70% cloudy, was <40%) and added a mosaic fallback to fill most such gaps automatically.');
  print('--- MODIS-Aqua (Chlorophyll-a, SST, Turbidity proxy) ---');
  print('Real NASA ocean-colour data, but coarse (~4km/pixel) - the fjord is narrower than that in many places, so these often show as one or two near-uniform colour blocks.');
  print('--- The "circle" ---');
  print('Every satellite layer is clipped to a 15km circle around GF3 (the white 1km ring is the same center, drawn smaller). Toggling a layer shows/hides that clip boundary - not a data artifact.');
  print('--- Honest note ---');
  print('If it is currently polar winter at 64N, a live fetch may return sparse or empty optical data (darkness/ice/cloud) - a real limitation of optical satellites this far north, not a bug.');
});
panel.add(satBtn);


// Forward-declared (assigned later, in the staged tail section below) so
// this button's closure can reference them once they exist.
var depthLayer, s2Layer, chlorALayer, sstLayer, kd490Layer, appMap;
var chlorA, sst, kd490, s2; // the ORIGINAL fixed Jun-Sep 2024 composites - kept
                            // around (not overwritten) so a Revert button can
                            // restore them after a live fetch swaps them out.

panel.add(lbl('WHAT THIS BUTTON DOES AND DOES NOT DO: it only refreshes the 4 satellite MAP LAYERS below (the pixels you see on the map) to a real 90-day window ending ~45 days ago (not literally today - satellite ocean-colour data has real processing latency; a real test showed asking for data ending today returns zero images). It does NOT change what the COUPLING ENGINE above uses - that is controlled separately by each variable\'s own In-situ/Satellite dropdown. Clicking this never makes the Coupling Engine "switch to satellite" by itself.', 8, '#5d4037'));
var liveFetchBtn = actionBtn('\uD83D\uDD04 Fetch LIVE Satellite Data (recent 90-day window) - MAP LAYERS ONLY', '#fff3cd', '#856404', '#ffc107', null);
panel.add(liveFetchBtn);
var liveFetchStatus = dynLbl('Layers start OFF for fast load. Enable via native Layers panel, or click above for a live fetch.', '#856404');
panel.add(liveFetchStatus);

liveFetchBtn.onClick(function() {
  liveFetchStatus.setValue('Fetching live composite from Earth Engine...');
  liveFetchStatus.style().set('color', '#f9a825');
  // FIX: a real test showed "last 90 days ending TODAY" returned ZERO
  // MODIS-Aqua images - not a bug, an expected real-world limitation.
  // Satellite Level-3 ocean-colour products have genuine processing
  // latency (days to weeks) before recent data is published to the
  // archive. Shifted the whole window back by 45 days so it targets
  // data that's actually had time to be processed, instead of asking
  // for a gap that hasn't been filled in yet.
  var PROCESSING_LATENCY_DAYS = 45;
  var today = new Date();
  var end = new Date(today.getTime() - PROCESSING_LATENCY_DAYS * 24 * 3600 * 1000);
  var endStr = end.toISOString().slice(0, 10);
  var past = new Date(end.getTime() - 90 * 24 * 3600 * 1000);
  var startStr = past.toISOString().slice(0, 10);
  print('=== Fetching LIVE satellite composite: ' + startStr + ' to ' + endStr + ' (shifted back ' + PROCESSING_LATENCY_DAYS + ' days from today to avoid the unprocessed-data gap) ===');

  var liveOceanColor = ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(startStr, endStr).filterBounds(REGION);
  var liveChlorA = liveOceanColor.select('chlor_a').mean().clip(REGION);
  var liveSst = liveOceanColor.select('sst').mean().clip(REGION);
  var liveKd490 = liveOceanColor.select('Kd_490').mean().clip(REGION);
  liveOceanColor.size().evaluate(function(count, err) {
    if (err) { print('[diagnostic] Could not check live MODIS-Aqua image count: ' + err); return; }
    print('[diagnostic] MODIS-Aqua L3SMI images available for ' + startStr + ' to ' + endStr + ' within 15km of GF3: ' + count);
    if (count === 0) print('[diagnostic] *** ZERO images in this live window too - confirms the region/dataset genuinely has no coverage here, not a one-off. ***');
  });
  var liveS2Collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterDate(startStr, endStr).filterBounds(REGION)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 70));
  var liveS2 = liveS2Collection.median().unmask(liveS2Collection.mosaic()).clip(REGION);

  chlorALayer.setEeObject(liveChlorA);
  sstLayer.setEeObject(liveSst);
  kd490Layer.setEeObject(liveKd490);
  s2Layer.setEeObject(liveS2);

  liveFetchStatus.setValue('Live composite loaded: ' + startStr + ' to ' + endStr + ' (fetched just now).');
  liveFetchStatus.style().set('color', '#1a4a2a');
});

// Revert button - closes a real gap: before this, once you clicked Fetch
// LIVE, there was no way back to the original fixed Jun-Sep 2024
// composite except re-running the whole script. This just calls
// setEeObject() with the ORIGINAL images (kept in the forward-declared
// chlorA/sst/kd490/s2 above, never overwritten - only the LAYER's
// displayed image was ever swapped, not these underlying variables).
var revertFetchBtn = actionBtn('\u21A9 Revert to original Jun-Sep 2024 composite', '#e0e0e0', '#333333', '#9e9e9e', function() {
  chlorALayer.setEeObject(chlorA);
  sstLayer.setEeObject(sst);
  kd490Layer.setEeObject(kd490);
  s2Layer.setEeObject(s2);
  liveFetchStatus.setValue('Reverted to the original fixed 2024-06-01 to 2024-09-30 composite.');
  liveFetchStatus.style().set('color', '#333333');
  print('=== Reverted satellite map layers to original Jun-Sep 2024 composite ===');
});
panel.add(revertFetchBtn);

// =================================================================
// 7. DATA SOURCES
// =================================================================
panel.add(sHead('DATA SOURCES', '#37474f'));
panel.add(lbl('In-situ: GEM MarineBasis Nuuk CTD (Mortensen et al. 2022). Real turbidity, fluorescence, salinity, temperature, density at Surface/Mid/Deep, 2005-2024. Satellite: NASA MODIS-Aqua L3SMI, ESA Copernicus Sentinel-2, NOAA ETOPO1.', 8, '#888888'));
var liveDataStatusLbl = dynLbl('In-situ data: hardcoded snapshot (2005-2024). Checking for a live monthly-refresh asset...', '#888888');
panel.add(liveDataStatusLbl);

// =================================================================
// 8-10. SATELLITE LAYERS + MAP + COMPOSE — genuinely staged this time.
// FIX: everything after the sidebar previously ran as ONE uninterrupted
// synchronous block. Browsers typically don't paint/flush console
// output until the current synchronous block yields back to the event
// loop - so if ui.root.add() (the final DOM mount of ~150+ widgets)
// was the real slow step, NONE of the print() checkpoints above it
// would ever become visible, even though they'd already executed.
// That exactly matches the report of a totally blank Console. Real
// fix: chain each remaining phase through ui.util.setTimeout(fn, 50),
// which genuinely returns control to the browser between phases, so
// each print() actually gets painted before the next (possibly slow)
// phase begins. This will tell us, for real this time, whether the
// mount step specifically is the bottleneck.
// =================================================================
print('[load 3/5] Sidebar built. Yielding to browser before defining satellite layer images...');

ui.util.setTimeout(function() {
  print('[load 3/5 running] Defining satellite layer images (lazy, server-side - no pixels computed yet)...');
  var recentStart = '2024-06-01';
  var recentEnd = '2024-09-30';

  var oceanColor = ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(recentStart, recentEnd).filterBounds(REGION);
  chlorA = oceanColor.select('chlor_a').mean().clip(REGION);
  sst = oceanColor.select('sst').mean().clip(REGION);
  kd490 = oceanColor.select('Kd_490').mean().clip(REGION);

  // REAL diagnostic, not a guess: query the actual number of MODIS-Aqua
  // images that pass the date+region filter for the FIXED composite.
  // If this comes back 0 (or very low), that directly confirms "the
  // source collection has no real coverage here for this window" as
  // the cause of Chlorophyll/SST/Turbidity never rendering - rather
  // than continuing to guess at visualization parameters or rendering
  // bugs. Runs once, async, doesn't block anything else loading.
  oceanColor.size().evaluate(function(count, err) {
    if (err) { print('[diagnostic] Could not check MODIS-Aqua image count: ' + err); return; }
    print('[diagnostic] MODIS-Aqua L3SMI images available for ' + recentStart + ' to ' + recentEnd + ' within 15km of GF3: ' + count);
    if (count === 0) {
      print('[diagnostic] *** ZERO images found. This is why Chlorophyll/SST/Turbidity never render - .mean() over an empty collection produces a fully transparent image. The fixed composite window needs to be widened or moved. ***');
    }
  });
  // Second, more precise check: even with nonzero SCENE count above,
  // every one of those scenes could still be 100% cloud/ice-masked
  // specifically over this small 15km circle, leaving zero VALID
  // pixels in the final .mean() composite despite images existing
  // nearby. This directly counts real unmasked pixels, which is the
  // actual determinant of whether anything renders.
  chlorA.reduceRegion({reducer: ee.Reducer.count(), geometry: REGION, scale: 4000, maxPixels: 1e9, bestEffort: true}).get('chlor_a').evaluate(function(pxCount, err2) {
    if (err2) { print('[diagnostic] Could not check valid pixel count: ' + err2); return; }
    print('[diagnostic] Valid (unmasked) chlor_a pixels within the study region: ' + pxCount);
    if (!pxCount || pxCount === 0) {
      print('[diagnostic] *** ZERO valid pixels even though scenes may exist nearby - every scene was cloud/ice-masked over this exact small area for this window. This is the direct cause of nothing rendering. ***');
    }
  });
  // Third check: real DATA EXISTS (confirmed above) but may still never
  // render if the actual values fall outside the 0-5 mg/m3 stretch the
  // legend/palette assumes. This is a real, documented failure mode:
  // standard ocean-colour chlorophyll algorithms (including MODIS-
  // Aqua's) are calibrated for open "Case 1" ocean water and are known
  // to misbehave in turbid, glacial-sediment-laden fjords like this one
  // - producing values that are negative, near-zero, or far above 5.
  // If real values sit entirely outside [0,5], the whole layer paints
  // as one flat color (usually the darkest end of the palette), which
  // can be visually indistinguishable from "nothing rendered" against
  // a dark ocean basemap.
  chlorA.reduceRegion({reducer: ee.Reducer.minMax(), geometry: REGION, scale: 4000, maxPixels: 1e9, bestEffort: true}).evaluate(function(stats, err3) {
    if (err3) { print('[diagnostic] Could not check chlor_a value range: ' + err3); return; }
    print('[diagnostic] Real chlor_a value range in this composite: min=' + stats.chlor_a_min + ', max=' + stats.chlor_a_max + ' (legend/palette currently assumes 0 to 5 mg/m3)');
    if (stats.chlor_a_min < -0.5 || stats.chlor_a_max > 20 || stats.chlor_a_max < 0.1) {
      print('[diagnostic] *** Real values fall well outside the assumed 0-5 range - this fjord\'s turbid, sediment-laden water is a documented case where standard ocean-colour chlorophyll algorithms misbehave. The palette stretch needs to match these real numbers, not a generic open-ocean assumption. ***');
    }
  });

  // Raised cloud threshold + mosaic fallback fixes a real artifact: a
  // strict median composite can leave a pixel with ZERO valid
  // observations (visible as a black/no-data ring, usually where
  // Sentinel-2 orbit swaths only partially overlap the study circle).
  // .unmask(mosaic) fills any such gap with the best single available
  // pixel instead of leaving it empty.
  var s2Collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterDate(recentStart, recentEnd).filterBounds(REGION)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 70));
  s2 = s2Collection.median().unmask(s2Collection.mosaic()).clip(REGION);

  var depthRaw = ee.Image('NOAA/NGDC/ETOPO1').select('bedrock').clip(REGION);
  var depthPositive = depthRaw.multiply(-1);
  var depthClass = ee.Image(0)
    .where(depthPositive.gte(0).and(depthPositive.lt(20)), 1)
    .where(depthPositive.gte(20).and(depthPositive.lt(300)), 2)
    .where(depthPositive.gte(300), 3)
    .updateMask(depthRaw.lt(0));
  // visDepthClass/visChl/visSst/visKd/visS2 are NOT redefined here -
  // they're the same objects declared early (near STUDY SITE), shared
  // by closure with the LEGEND section built earlier in the sidebar.
  // This guarantees the legend and these map layers can never disagree.

  print('[load 4/5] Satellite layer images defined. Yielding before building the Map + registering layers...');

  ui.util.setTimeout(function() {
    print('[load 4/5 running] Building map, registering layers (all default OFF)...');
    appMap = ui.Map();
    appMap.setOptions('SATELLITE');
    appMap.centerObject(NUUK, 10);
    appMap.style().set({cursor: 'crosshair'});

    // FIX: a real error - "Layer error: Map must be specified" - showed
    // up on every layer the Fetch LIVE Satellite Data button tries to
    // update. Root cause: Map.addLayer()'s shorthand return value does
    // NOT reliably stay bound to this Map for a later .setEeObject()
    // call (my earlier mock test harness assumed it would, which is
    // exactly the kind of thing a mock can't catch - it only proves the
    // code I WROTE behaves as I assumed, not that the assumption
    // matches real Earth Engine behavior). The documented, reliable
    // pattern is to construct the ui.Map.Layer object directly first,
    // THEN register it via layers().add() - only layers built this way
    // can safely have setEeObject() called on them later.
    depthLayer = ui.Map.Layer(depthClass, visDepthClass, 'Depth (Surface/Mid/Deep) [satellite: ETOPO1, coarse]', false, 0.65);
    s2Layer = ui.Map.Layer(s2, visS2, 'Sentinel-2 true colour photo [satellite: Sentinel-2]', false, 1);
    chlorALayer = ui.Map.Layer(chlorA, visChl, 'Chlorophyll-a [satellite: MODIS-Aqua, coarse]', false, 0.75);
    sstLayer = ui.Map.Layer(sst, visSst, 'SST, deg C [satellite: MODIS-Aqua, coarse]', false, 0.75);
    kd490Layer = ui.Map.Layer(kd490, visKd, 'Turbidity proxy Kd_490 [satellite: MODIS-Aqua, coarse]', false, 0.75);
    appMap.layers().add(depthLayer);
    appMap.layers().add(s2Layer);
    appMap.layers().add(chlorALayer);
    appMap.layers().add(sstLayer);
    appMap.layers().add(kd490Layer);
    appMap.addLayer(ringFC.style({color: 'FFFFFF', fillColor: '00000000', width: 2}), {}, 'Study zone [in-situ] (1km ring)');
    appMap.addLayer(siteFC.style({color: '2e7d32', pointSize: 8}), {}, 'GF3 station [in-situ], Godthaabsfjord');

    var mapInfoBox = ui.Panel({
      style: {position: 'top-left', padding: '8px', backgroundColor: 'rgba(6, 20, 12, 0.92)', border: '1px solid #2e7d32', width: '220px'}
    });
    mapInfoBox.add(ui.Label('SW Greenland / Nuuk', {color: '#66bb6a', backgroundColor: 'rgba(0,0,0,0)', fontWeight: 'bold', fontSize: '13px', margin: '0 0 2px 0'}));
    mapInfoBox.add(ui.Label('IN-SITU CTD | GF3, real depth to 392 dbar', {color: '#81c784', backgroundColor: 'rgba(0,0,0,0)', fontSize: '10px'}));
    mapInfoBox.add(ui.Label('64.117 N, -51.883 W', {color: '#ffffff', backgroundColor: 'rgba(0,0,0,0)', fontWeight: 'bold', fontSize: '10px', margin: '2px 0 0 0'}));
    var etopoLabel = ui.Label('ETOPO1 depth here: querying live...', {color: '#ffeb3b', backgroundColor: 'rgba(0,0,0,0)', fontSize: '9px', margin: '4px 0 0 0'});
    mapInfoBox.add(etopoLabel);
    appMap.add(mapInfoBox);

    // Real, single, cheap Earth Engine query - the actual ETOPO1 value
    // AT the exact GF3 coordinates, not an eyeballed color. Answers
    // "is the map really showing this point as shallow?" with a number.
    depthRaw.reduceRegion({reducer: ee.Reducer.first(), geometry: NUUK, scale: 500, maxPixels: 1e6}).get('bedrock').evaluate(function(val, err) {
      if (err || val === null || val === undefined) {
        etopoLabel.setValue('ETOPO1 depth here: query failed or no data at this exact point.');
        return;
      }
      var etopoDepthM = Math.round(-val);
      etopoLabel.setValue('ETOPO1 raw value here: ' + etopoDepthM + 'm (vs real CTD: 392 dbar). ' +
        (etopoDepthM < 300 ? 'Confirms the coarse-resolution mismatch - this exact 1.85km pixel reads shallower than the real station.' : 'Roughly agrees with the real CTD depth at this point.'));
    });

    print('[load 5/5] Map + layers built. Yielding before final DOM mount - THIS is the single most suspected slow step (one-time layout of the whole widget tree). If the Console stops here, this confirms it.');

    ui.util.setTimeout(function() {
      print('[load 5/5 running] Mounting to root now...');
      ui.root.clear();
      var rootSplit = ui.SplitPanel({firstPanel: panel, secondPanel: appMap, orientation: 'horizontal', wipe: false, style: {stretch: 'both'}});
      ui.root.add(rootSplit);
      print('=== MOUNT COMPLETE. NUUK TOOL v5 FULLY LOADED AND VISIBLE. Any hang after this line is inside a specific button click, not page load. ===');

      // Live-data check runs AFTER the page is already up and visible -
      // deliberately not blocking page load. If the pipeline's asset
      // doesn't exist yet, this is a quick, cheap failure; if it does,
      // NUUK_DATA gets updated in place and any Coupling Engine run
      // after this point uses the refreshed values automatically.
      tryLoadLiveNuukData(function(res) {
        if (res.loaded) {
          liveDataStatusLbl.setValue('In-situ data: LIVE asset merged (' + res.monthsInAsset + ' months from the monthly pipeline, ' + res.seriesTouched + ' series updated).');
          liveDataStatusLbl.style().set('color', '#1a4a2a');
          print('[live-data] Merged ' + res.monthsInAsset + ' months from ' + LIVE_ASSET_ID + ' into NUUK_DATA.');
        } else {
          liveDataStatusLbl.setValue('In-situ data: hardcoded snapshot (2005-2024). ' + res.reason);
          liveDataStatusLbl.style().set('color', '#888888');
          print('[live-data] No live asset used: ' + res.reason);
        }
      });
    }, 50);
  }, 50);
}, 50);
