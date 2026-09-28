/* =============================================================
 * MolBench · core.js
 * 分子生物学核心算法库（零依赖 / 纯前端 / 可离线）
 * -------------------------------------------------------------
 * 覆盖范围：
 *   1. FASTA 解析与序列类型识别
 *   2. 碱基组成 / GC 含量 / GC 滑动窗口
 *   3. 反向互补 / 转录 / 翻译 / 六框 ORF 查找
 *   4. Tm：Wallace 规则、GC% 经验式、最近邻法(SantaLucia 1998)
 *   5. 分子量（ssDNA/ssRNA/dsDNA/蛋白）与 OD260 换算
 *   6. 限制性内切酶位点查找、线性/环状酶切、双酶切片段
 *   7. 引物评估：GC clamp、发夹、自身/交叉二聚体、PCR 产物
 * ============================================================= */
(function (global) {
  'use strict';

  var WATER = 18.0153;

  /* ---------- 标准遗传密码表 (NCBI translation table 1) ---------- */
  var CODON_TABLE = {
    TTT: 'F', TTC: 'F', TTA: 'L', TTG: 'L',
    CTT: 'L', CTC: 'L', CTA: 'L', CTG: 'L',
    ATT: 'I', ATC: 'I', ATA: 'I', ATG: 'M',
    GTT: 'V', GTC: 'V', GTA: 'V', GTG: 'V',
    TCT: 'S', TCC: 'S', TCA: 'S', TCG: 'S',
    CCT: 'P', CCC: 'P', CCA: 'P', CCG: 'P',
    ACT: 'T', ACC: 'T', ACA: 'T', ACG: 'T',
    GCT: 'A', GCC: 'A', GCA: 'A', GCG: 'A',
    TAT: 'Y', TAC: 'Y', TAA: '*', TAG: '*',
    CAT: 'H', CAC: 'H', CAA: 'Q', CAG: 'Q',
    AAT: 'N', AAC: 'N', AAA: 'K', AAG: 'K',
    GAT: 'D', GAC: 'D', GAA: 'E', GAG: 'E',
    TGT: 'C', TGC: 'C', TGA: '*', TGG: 'W',
    CGT: 'R', CGC: 'R', CGA: 'R', CGG: 'R',
    AGT: 'S', AGC: 'S', AGA: 'R', AGG: 'R',
    GGT: 'G', GGC: 'G', GGA: 'G', GGG: 'G'
  };

  /* ---------- 氨基酸残基分子量 (Da, 链内残基, 已脱水) ---------- */
  var AA_RESIDUE = {
    A: 71.0788, R: 156.1875, N: 114.1038, D: 115.0886, C: 103.1388,
    E: 129.1155, Q: 128.1307, G: 57.0519, H: 137.1411, I: 113.1594,
    L: 113.1594, K: 128.1741, M: 131.1926, F: 147.1766, P: 97.1167,
    S: 87.0782, T: 101.1051, W: 186.2132, Y: 163.1760, V: 99.1326
  };

  /* ---------- 核苷酸残基分子量 (Da, 链内残基, 已脱水) ---------- */
  var DNA_RESIDUE = { A: 313.21, C: 289.18, G: 329.21, T: 304.20 };
  var RNA_RESIDUE = { A: 329.21, C: 305.18, G: 345.21, U: 306.17 };

  /* ---------- 简并碱基 ---------- */
  var IUPAC = {
    A: 'A', C: 'C', G: 'G', T: 'T', U: 'U',
    R: '[AG]', Y: '[CT]', S: '[GC]', W: '[AT]', K: '[GT]', M: '[AC]',
    B: '[CGT]', D: '[AGT]', H: '[ACT]', V: '[ACG]', N: '[ACGTU]'
  };

  /* ---------- 互补碱基（含简并） ---------- */
  var COMPLEMENT = {
    A: 'T', T: 'A', G: 'C', C: 'G', U: 'A',
    R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K',
    B: 'V', V: 'B', D: 'H', H: 'D', N: 'N'
  };

  /* ---------- 最近邻热力学参数 (SantaLucia 1998, unified) ----------
   * key = 5'->3' 相邻二核苷酸; value = [dH (kcal/mol), dS (cal/mol·K)]
   * 反向互补的二核苷酸共用一组参数，故 16 种组合全部覆盖。
   * ---------------------------------------------------------------- */
  var NN = {
    AA: [-7.9, -22.2], TT: [-7.9, -22.2],
    AT: [-7.2, -20.4], TA: [-7.2, -21.3],
    CA: [-8.5, -22.7], TG: [-8.5, -22.7],
    GT: [-8.4, -22.4], AC: [-8.4, -22.4],
    CT: [-7.8, -21.0], AG: [-7.8, -21.0],
    GA: [-8.2, -22.2], TC: [-8.2, -22.2],
    CG: [-10.6, -27.2], GC: [-9.8, -24.4],
    GG: [-8.0, -19.9], CC: [-8.0, -19.9]
  };
  var NN_INIT_GC = [0.1, -2.8];
  var NN_INIT_AT = [2.3, 4.1];
  var R_GAS = 1.9872; // cal/(mol·K)

  /* =========================================================
   * 1. 解析与识别
   * ========================================================= */
  /** 识别输入格式，供 parseFasta 分派 */
  function detectFormat(text) {
    if (/^\s*>/.test(text)) return 'fasta';
    if (/^LOCUS\s+\S/m.test(text) || /^\s*ORIGIN\s*$/m.test(text)) return 'genbank';
    if (/^ID\s+\S/m.test(text) || /^SQ\s/m.test(text)) return 'embl';
    return 'raw';
  }

  /**
   * 解析 GenBank 格式：只取 ORIGIN 段之后的序列。
   * 若不做此处理，LOCUS / DEFINITION / FEATURES 里的英文单词会被当成碱基，
   * 产生完全错误的"序列"（实测 40 nt 的记录会被解析成 158 nt 垃圾）。
   */
  function parseGenBank(text) {
    var lines = String(text).split(/\r?\n/);
    var name = '', seq = '', inOrigin = false, i, line;
    for (i = 0; i < lines.length; i++) {
      line = lines[i];
      if (/^LOCUS\s+/.test(line)) {
        var parts = line.trim().split(/\s+/);
        if (parts[1]) name = parts[1];
      } else if (/^DEFINITION\s+/.test(line)) {
        var def = line.replace(/^DEFINITION\s+/, '').trim();
        if (def && !name) name = def;
        else if (def && name) name = name; // LOCUS 名优先
      } else if (/^\s*ORIGIN/.test(line)) {
        inOrigin = true;
      } else if (/^\s*\/\//.test(line)) {
        inOrigin = false;
      } else if (inOrigin) {
        seq += line.replace(/[^A-Za-z]/g, '');
      }
    }
    seq = seq.toUpperCase();
    return { name: name, seq: seq, type: detectType(seq), format: 'genbank', features: parseGenBankFeatures(text) };
  }

  /**
   * 解析 GenBank 的 FEATURES 段，提取 CDS / gene / mRNA 等注释。
   * 支持 complement(...)、join(...)、<1 / >100 等位置写法，以及跨多行的 qualifier
   * （尤其是 /translation 常常折成多行）。
   */
  function parseGenBankFeatures(text) {
    var lines = String(text).split(/\r?\n/);
    var feats = [];
    var inFeat = false, cur = null, curKey = null;
    var i, line, m, q;

    for (i = 0; i < lines.length; i++) {
      line = lines[i];
      if (/^\s*FEATURES/.test(line)) { inFeat = true; continue; }
      if (!inFeat) continue;
      if (/^\s*(ORIGIN|CONTIG|BASE COUNT)\b/.test(line)) break;
      if (/^\s*\/\//.test(line)) break;
      if (!line.trim()) continue;

      // 新的 feature 行：缩进 5 空格 + 类型 + 位置
      m = line.match(/^ {1,10}([A-Za-z_][A-Za-z0-9_'-]*)\s+(\S.*)$/);
      if (m && m[1].indexOf('/') !== 0) {
        if (cur) feats.push(cur);
        cur = { type: m[1], location: m[2].trim(), qualifiers: {} };
        curKey = null;
        continue;
      }
      // qualifier 行
      q = line.match(/^\s*\/([A-Za-z_][A-Za-z0-9_]*)\s*=?\s*(.*)$/);
      if (q && cur) {
        var key = q[1], rest = q[2];
        if (rest.charAt(0) === '"') rest = rest.slice(1);
        var closed = /"$/.test(rest);
        if (closed) rest = rest.slice(0, -1);
        cur.qualifiers[key] = rest;
        curKey = closed ? null : key;   // 引号未闭合 → 后续行继续追加
        continue;
      }
      // qualifier 的续行
      if (cur && curKey) {
        var add = line.trim();
        if (/"$/.test(add)) { add = add.slice(0, -1); cur.qualifiers[curKey] += add; curKey = null; }
        else cur.qualifiers[curKey] += add;
      }
    }
    if (cur) feats.push(cur);

    // 位置解析
    feats.forEach(function (f) {
      var loc = parseLocation(f.location);
      f.strand = loc.strand;
      f.segments = loc.segments;
      f.join = loc.join;
      if (loc.segments.length) {
        f.start = loc.segments[0].start;
        f.end = loc.segments[loc.segments.length - 1].end;
        // join 的总长是各段之和，不是首尾跨度（否则内含子会被算进去）
        f.length = loc.segments.reduce(function (t, sg) {
          return t + Math.abs(sg.end - sg.start) + 1;
        }, 0);
      }
    });
    return feats;
  }

  /**
   * 解析 GenBank location 字符串。
   * 支持 `123..456`、单碱基 `123`、`complement(...)`、`join(a..b, c..d)`、
   * 以及 `<1..>100` 这类不确定端点（尖括号被忽略）。
   * 之前的实现只取首尾数字，`join(1..30,50..100)` 会被当成 1..100，
   * 把内含子也算进 CDS 长度里。
   */
  function parseLocation(location) {
    var raw = String(location || '');
    var strand = /complement/.test(raw) ? '-' : '+';
    // 去掉 `<` / `>`（GenBank 用它表示端点不确定），否则 `<1..>100` 会被拆成
    // "1" 和 "100" 两个单碱基段，而不是一个 1..100 的区间
    var loc = raw.replace(/[<>]/g, '');
    var segments = [];
    var re = /(\d+)\s*\.\.\s*(\d+)|(\d+)/g;
    var m;
    while ((m = re.exec(loc)) !== null) {
      if (m[1] !== undefined) segments.push({ start: parseInt(m[1], 10), end: parseInt(m[2], 10) });
      else segments.push({ start: parseInt(m[3], 10), end: parseInt(m[3], 10) });
    }
    return { strand: strand, segments: segments, join: /join|order/.test(loc) };
  }

  /**
   * 按 feature 位置提取子序列（自动处理负链与 join 多段）。
   * join 的处理顺序：各段按坐标升序拼接，再整体反向互补（若为负链）——
   * 这正是 GenBank 对 `complement(join(...))` 的语义。
   */
  function extractFeature(seq, feat) {
    if (!seq || !feat) return '';
    var s = String(seq).toUpperCase();
    var segs = (feat.segments && feat.segments.length)
      ? feat.segments
      : (feat.start && feat.end ? [{ start: feat.start, end: feat.end }] : []);
    if (!segs.length) return '';
    var parts = [];
    segs.forEach(function (sg) {
      var a = Math.max(0, Math.min(sg.start, sg.end) - 1);
      var b = Math.min(s.length, Math.max(sg.start, sg.end));
      if (b > a) parts.push(s.slice(a, b));
    });
    var joined = parts.join('');
    return feat.strand === '-' ? revComp(joined) : joined;
  }

  /** 解析 EMBL 格式：取 SQ 段到 // 之间的序列 */
  function parseEmbl(text) {
    var lines = String(text).split(/\r?\n/);
    var name = '', seq = '', inSq = false, i, line;
    for (i = 0; i < lines.length; i++) {
      line = lines[i];
      if (/^ID\s+/.test(line)) {
        var first = line.replace(/^ID\s+/, '').trim().split(/[;\s]/)[0];
        if (first) name = first;
      } else if (/^DE\s+/.test(line)) {
        var de = line.replace(/^DE\s+/, '').trim();
        if (de && !name) name = de;
      } else if (/^SQ\s/.test(line)) {
        inSq = true;
      } else if (/^\s*\/\//.test(line)) {
        inSq = false;
      } else if (inSq) {
        seq += line.replace(/[^A-Za-z]/g, '');
      }
    }
    seq = seq.toUpperCase();
    return { name: name, seq: seq, type: detectType(seq), format: 'embl' };
  }

  function parseFasta(raw) {
    if (!raw) return { name: '', seq: '', type: 'unknown' };
    var text = String(raw);
    var fmt = detectFormat(text);
    // 带 ORIGIN / SQ 段的结构化格式优先专用解析器；解析不到序列时再回退到通用清洗
    if (fmt === 'genbank') {
      var gb = parseGenBank(text);
      if (gb.seq) return gb;
    } else if (fmt === 'embl') {
      var em = parseEmbl(text);
      if (em.seq) return em;
    }

    var name = '';
    var body = text;

    if (text.charAt(0) === '>') {
      var nl = text.indexOf('\n');
      if (nl === -1) { name = text.slice(1).trim(); body = ''; }
      else { name = text.slice(1, nl).trim(); body = text.slice(nl + 1); }
    }
    // 去掉空白、数字、行号、常见比对符号
    var seq = body.replace(/\s+/g, '').replace(/[0-9]/g, '').replace(/[-*.]/g, '').toUpperCase();

    return { name: name, seq: seq, type: detectType(seq), format: 'fasta' };
  }

  function detectType(seq) {
    if (!seq) return 'unknown';
    var s = seq.toUpperCase();
    var n = s.length;
    var nt = 0, aa = 0;
    var aaOnly = 'DEFHIKLMNPQRSVWY'; // 绝不会出现在核酸里的残基字母
    for (var i = 0; i < n; i++) {
      var c = s.charAt(i);
      if ('ACGTUNRYKMSWBDHV'.indexOf(c) !== -1) nt++;
      else if (aaOnly.indexOf(c) !== -1) aa++;
    }
    if (aa > n * 0.02) return 'protein';
    if (nt === 0) return 'unknown';
    // 有 U 且几乎无 T → RNA
    var u = (s.match(/U/g) || []).length;
    var t = (s.match(/T/g) || []).length;
    if (u > 0 && t === 0) return 'rna';
    return 'dna';
  }

  function isValidNucleic(seq) {
    if (!seq) return false;
    var s = String(seq).toUpperCase();
    for (var i = 0; i < s.length; i++) {
      if (!IUPAC[s.charAt(i)]) return false;
    }
    return true;
  }

  /* =========================================================
   * 2. 组成与 GC
   * ========================================================= */
  function baseCount(seq) {
    var s = String(seq).toUpperCase();
    var c = { A: 0, C: 0, G: 0, T: 0, U: 0, N: 0, other: 0 };
    for (var i = 0; i < s.length; i++) {
      var ch = s.charAt(i);
      if (c[ch] === undefined) c.other++;
      else c[ch]++;
    }
    return c;
  }

  function gcContent(seq) {
    if (!seq) return 0;
    var c = baseCount(seq);
    var gc = c.G + c.C;
    var at = c.A + c.T + c.U;
    var tot = gc + at;
    return tot === 0 ? 0 : (gc / tot) * 100;
  }

  function isGC(c) { return c === 'G' || c === 'C'; }
  function isAT(c) { return c === 'A' || c === 'T' || c === 'U'; }

  /**
   * GC 滑动窗口：返回 [{pos, gc}]，pos 为窗口起点(1-based)
   * 用增量滑动把复杂度从 O(n×win) 降到 O(n)：每次只处理移出/移入的 step 个字符。
   * 分母只统计 GC+AT（与 gcContent 语义一致，N 等简并字符不计入）。
   */
  function gcWindows(seq, win, step) {
    win = win || 50; step = step || 10;
    var src = String(seq).toUpperCase();
    var out = [];
    var n = src.length;
    if (n < win || win < 1) return out;

    var gc = 0, valid = 0, i;
    for (i = 0; i < win; i++) {
      var c = src.charAt(i);
      if (isGC(c)) { gc++; valid++; }
      else if (isAT(c)) valid++;
    }
    out.push({ pos: 1, gc: valid ? (gc / valid) * 100 : 0 });

    for (var start = step; start + win <= n; start += step) {
      for (var k = 0; k < step; k++) {
        var oc = src.charAt(start - step + k);
        if (isGC(oc)) { gc--; valid--; } else if (isAT(oc)) valid--;
        var ic = src.charAt(start + win - step + k);
        if (isGC(ic)) { gc++; valid++; } else if (isAT(ic)) valid++;
      }
      out.push({ pos: start + 1, gc: valid ? (gc / valid) * 100 : 0 });
    }
    return out;
  }

  /* =========================================================
   * 3. 反向互补 / 转录 / 翻译 / ORF
   * ========================================================= */
  function revComp(seq) {
    var s = String(seq).toUpperCase();
    var out = '';
    for (var i = s.length - 1; i >= 0; i--) {
      out += COMPLEMENT[s.charAt(i)] || 'N';
    }
    return out;
  }

  function transcribe(dna) { return dna.replace(/T/g, 'U'); }
  function reverseTranscribe(rna) { return rna.replace(/U/g, 'T'); }

  function translate(ntSeq, frameOffset) {
    frameOffset = frameOffset || 0;
    // 核心 API 必须自己对大小写设防：直接查密码表时小写会全部落到 X
    var s = String(ntSeq).toUpperCase().replace(/U/g, 'T');
    var out = '';
    for (var i = frameOffset; i + 2 < s.length; i += 3) {
      var codon = s.substr(i, 3);
      out += CODON_TABLE[codon] || 'X';
    }
    return out;
  }

  /**
   * 构造一条 ORF 记录。
   * 两个坐标必须分清：
   *  - chainStart/chainEnd：在该链自身坐标系里的位置（负链是 revComp 上的位置），
   *    环状拓扑做环形映射时需要它；
   *  - ntStart/ntEnd：**统一映射回原序列**的位置。
   * 负链若不映射，用户拿这个区间去原序列取子序列会取到完全错误的位置。
   */
  function makeOrf(st, frameIdx, aaStart, aaStop, prot, n, proteinOverride) {
    var chainStart = frameIdx + aaStart * 3 + 1;      // 1-based，含终止密码子
    var chainEnd = frameIdx + (aaStop + 1) * 3;
    var isMinus = st.strand === '-';
    return {
      strand: st.strand,
      frame: frameIdx + 1,
      chainStart: chainStart,
      chainEnd: chainEnd,
      ntStart: isMinus ? n - chainEnd + 1 : chainStart,
      ntEnd: isMinus ? n - chainStart + 1 : chainEnd,
      aaLength: aaStop - aaStart,
      protein: proteinOverride || prot.slice(aaStart, aaStop),
      // 核酸长度统一为 (aa + 1) × 3，即包含终止密码子
      nt: st.seq.substr(frameIdx + aaStart * 3, (aaStop - aaStart + 1) * 3)
    };
  }

  /**
   * 六框 ORF 查找
   * @param {string} seq 核酸序列（DNA，可含 U）
   * @param {number} minAa 最小氨基酸长度（不含终止符）
   * @param {boolean} requireStart 是否要求 ATG 起始
   * @returns {Array} ORF 列表，按长度降序
   */
  function findORFs(seq, minAa, requireStart) {
    // 不能写成 minAa || 50：那样用户传入 0 会被静默替换成默认值 50。
    // 同时把下限收到 1，避免 minAa=0 时在长序列上产出海量记录。
    minAa = (minAa === undefined || minAa === null || isNaN(minAa)) ? 50 : Math.max(1, minAa);
    requireStart = requireStart !== false;
    var dna = seq.replace(/U/g, 'T').toUpperCase();
    var n = dna.length;
    var rc = revComp(dna);
    var out = [];
    var strands = [
      { seq: dna, strand: '+', label: '正链' },
      { seq: rc, strand: '-', label: '负链' }
    ];

    for (var s = 0; s < strands.length; s++) {
      var st = strands[s];
      for (var f = 0; f < 3; f++) {
        var prot = translate(st.seq, f);
        var i = 0;
        if (requireStart) {
          // 单遍扫描，维护"尚未闭合的起始密码子"。
          // 早期实现是"取第一个 M、跳到它后面的第一个 *"，遇到 M1…M2…* 这种嵌套时
          // 只会报出 M1 起始的长 ORF，而漏掉 M2 起始的短 ORF。
          // 现在每遇到一个终止子，就结算所有与它能构成合法长度的起始位点。
          var open = [];
          for (var p = 0; p < prot.length; p++) {
            var ch = prot.charAt(p);
            if (ch === 'M') open.push(p);
            else if (ch === '*') {
              for (var q = 0; q < open.length; q++) {
                if (p - open[q] < minAa) continue;   // 长度不足 minAa 的起始丢弃
                out.push(makeOrf(st, f, open[q], p, prot, n));
              }
              open.length = 0;   // 终止子之后，之前未闭合的起始都不再有效
            }
          }
        } else {
          while (i < prot.length) {
            // 不要求起始密码子：取两个终止子之间的最长片段
            var a = (i === 0) ? -1 : i - 1;
            var b = prot.indexOf('*', i);
            var seg = (b === -1) ? prot.slice(a + 1) : prot.slice(a + 1, b);
            if (seg.length >= minAa) {
              var segEnd = (b === -1 ? prot.length : b);
              out.push(makeOrf(st, f, a + 1, segEnd, prot, n, seg));
            }
            if (b === -1) break;
            i = b + 1;
          }
        }
      }
    }
    out.sort(function (a, b) { return b.aaLength - a.aaLength; });
    return out;
  }

  /**
   * 环状序列的六框 ORF 查找。
   * 环状 DNA 上的基因（或用户截取的片段）可能跨越首尾接缝，
   * 直接在线性序列上扫描会漏掉这类 ORF。
   * 做法：在 seq+seq 上扫描，只保留起点落在前半圈、长度不超过一圈的结果，
   * 再把坐标取模映射回原序列，并用 wrap 标记跨接缝。
   */
  function findORFsCircular(seq, minAa, requireStart) {
    var s = String(seq || '').replace(/U/g, 'T').toUpperCase();
    var n = s.length;
    if (!n) return [];
    if (n < 9) return findORFs(s, minAa, requireStart);   // 太短，环形绕不出一圈

    var raw = findORFs(s + s, minAa, requireStart);
    var N2 = 2 * n;
    var out = [], seen = {};

    raw.forEach(function (o) {
      var cs = o.chainStart, ce = o.chainEnd;
      if (cs > n) return;                 // 起点必须落在前半圈，避免重复计数
      if (ce - cs + 1 > n) return;        // ORF 不能超过一整圈
      var start, end, wrap;
      if (o.strand === '+') {
        start = cs;
        end = ((ce - 1) % n) + 1;
        wrap = ce > n;
      } else {
        // 负链坐标在 revComp(seq+seq) 上；先映射回 seq+seq，再取模回原序列
        var dStart = N2 - ce + 1;
        var dEnd = N2 - cs + 1;
        start = ((dStart - 1) % n) + 1;
        end = ((dEnd - 1) % n) + 1;
        wrap = dEnd > n;
      }
      var key = o.strand + ':' + start + ':' + end;
      if (seen[key]) return;
      seen[key] = 1;
      out.push({
        strand: o.strand,
        frame: o.frame,
        ntStart: start,
        ntEnd: end,
        aaLength: o.aaLength,
        protein: o.protein,
        nt: o.nt,
        wrap: wrap
      });
    });

    out.sort(function (a, b) { return b.aaLength - a.aaLength; });
    return out;
  }

  /* =========================================================
   * 4. 解链温度 Tm
   * ========================================================= */
  function tmWallace(seq) {
    var s = seq.replace(/U/g, 'T').toUpperCase();
    var c = baseCount(s);
    var at = c.A + c.T, gc = c.G + c.C;
    return 2 * at + 4 * gc;
  }

  function tmGC(seq) {
    var s = seq.replace(/U/g, 'T').toUpperCase();
    var n = s.length;
    if (n === 0) return 0;
    var c = baseCount(s);
    var gc = c.G + c.C;
    return 64.9 + 41 * (gc - 16.4) / n;
  }

  /**
   * 最近邻累加：返回完美配对双链的 [ΔH(kcal/mol), ΔS(cal/mol·K)]，
   * 已含两端 initiation 与盐校正。Tm 与 ΔG 共用这一份参数，避免两处各写一套。
   */
  function nnSum(seq, opts) {
    opts = opts || {};
    var s = seq.replace(/U/g, 'T').toUpperCase().replace(/[^ACGT]/g, '');
    if (s.length < 2) return null;
    var dH = 0, dS = 0;
    for (var i = 0; i < s.length - 1; i++) {
      var p = NN[s.substr(i, 2)];
      if (!p) return null;
      dH += p[0];
      dS += p[1];
    }
    // 两端 initiation
    var init5 = (s.charAt(0) === 'G' || s.charAt(0) === 'C') ? NN_INIT_GC : NN_INIT_AT;
    var init3 = (s.charAt(s.length - 1) === 'G' || s.charAt(s.length - 1) === 'C') ? NN_INIT_GC : NN_INIT_AT;
    dH += init5[0] + init3[0];
    dS += init5[1] + init3[1];

    // 盐校正（Owczarzy 2004 简化式）：把 Mg2+ 折算成等效单价阳离子
    // 注意：下列浓度均以 mM 为单位参与运算，最后统一折算为 M
    var naMM = (opts.na === undefined ? 50 : opts.na);
    var mgMM = (opts.mg === undefined ? 1.5 : opts.mg);
    var dntpMM = (opts.dntp === undefined ? 0.2 : opts.dntp);
    if (dntpMM > 0 && mgMM > 0) {
      // dNTP 螯合等量的 Mg2+，仅游离 Mg2+ 参与校正
      var freeMg = mgMM - dntpMM;
      mgMM = freeMg > 0.01 ? freeMg : 0.01;
    }
    var naEq = mgMM > 0 ? (naMM + 120 * Math.sqrt(mgMM)) / 1000 : naMM / 1000;
    if (!(naEq > 0)) naEq = 0.05;
    dS = dS + 0.368 * (s.length - 1) * Math.log(naEq);

    return [dH, dS];
  }

  /**
   * 最近邻法 Tm (SantaLucia 1998)
   * @param {string} seq
   * @param {object} opts {conc: 引物总浓度 nM, na: Na+ 浓度 mM, mg: Mg2+ mM, dntp: mM}
   */
  function tmNearestNeighbor(seq, opts) {
    opts = opts || {};
    var sum = nnSum(seq, opts);
    if (!sum) return NaN;
    // 引物总浓度 (M)。非自身互补且引物过量 → 用 CT/4
    var conc = (opts.conc === undefined ? 500 : opts.conc) * 1e-9; // nM → M
    var ct = conc / 4;
    if (!(ct > 0)) ct = 1.25e-7;
    return (sum[0] * 1000) / (sum[1] + R_GAS * Math.log(ct)) - 273.15;
  }

  /* =========================================================
   * 4b. 二级结构热力学（ΔG, 37 °C）
   * ========================================================= */
  var T37 = 310.15;

  /** ΔG(37 °C) = ΔH − T·ΔS/1000，单位 kcal/mol */
  function nnDeltaG37(seq, opts) {
    var sum = nnSum(seq, opts);
    if (!sum) return NaN;
    return sum[0] - T37 * sum[1] / 1000;
  }

  // 发夹环 initiation 自由能（kcal/mol, 37 °C），来自常用 DNA 折叠参数
  var HAIRPIN_LOOP = { 3: 5.4, 4: 5.6, 5: 5.7, 6: 5.4, 7: 6.0, 8: 5.5, 9: 6.4, 10: 6.5 };
  /** 环长罚值：≤10 查表，更长的按 Jacobson-Stockmayer 外推 */
  function loopPenalty(len) {
    if (HAIRPIN_LOOP[len] !== undefined) return HAIRPIN_LOOP[len];
    if (len < 3) return 8.0;
    return HAIRPIN_LOOP[10] + 1.079 * Math.log(len / 10);
  }

  /**
   * 发夹 ΔG：枚举环起止，求茎区最长连续完美配对，取最稳定结构。
   * 茎区用最近邻参数求和，加上环 initiation 罚值。
   * 注意：这是简化模型（不处理凸起/内部环/错配），用于引物快速筛查。
   * @returns {{dG:number, stem:number, loop:number, pos:number}}
   */
  function hairpinDeltaG(seq, opts) {
    var s = seq.replace(/U/g, 'T').toUpperCase().replace(/[^ACGT]/g, '');
    var n = s.length;
    var best = { dG: 0, stem: 0, loop: 0, pos: -1 };
    if (n < 9) return best;

    for (var i = 3; i < n - 3; i++) {
      for (var loop = 3; loop <= 30; loop++) {
        var j = i + loop;
        if (j >= n) break;
        var arm = '';
        for (var k = 0; k < 12; k++) {   // 茎长上限 12 bp：引物发夹不会更长，同时控制枚举开销
          var left = i - 1 - k, right = j + k;
          if (left < 0 || right >= n) break;
          if ((COMPLEMENT[s.charAt(left)] || 'N') !== s.charAt(right)) break;
          arm = s.charAt(left) + arm;   // 5′ 臂按 5′→3′ 方向拼装
        }
        if (arm.length < 3) continue;
        var stemG = arm.length >= 2 ? nnDeltaG37(arm, opts) : 1.5;
        if (isNaN(stemG)) continue;
        var dG = stemG + loopPenalty(loop);
        if (dG < best.dG) best = { dG: dG, stem: arm.length, loop: loop, pos: i };
      }
    }
    return best;
  }

  /**
   * 二聚体 ΔG：把 b 反向互补后与 a 平行比对，枚举所有偏移，
   * 对每段连续完美配对区（≥3）计算 ΔG，取最稳定者。
   * 两个引物互查即为交叉二聚体；a === b 即自身二聚体。
   * @returns {{dG:number, run:number, offset:number}}
   */
  function dimerDeltaG(a, b, opts) {
    var A = String(a || '').replace(/U/g, 'T').toUpperCase().replace(/[^ACGT]/g, '');
    var B = revComp(String(b || '').replace(/U/g, 'T').toUpperCase().replace(/[^ACGT]/g, ''));
    var best = { dG: 0, run: 0, offset: 0 };
    if (!A.length || !B.length) return best;

    for (var off = -(B.length - 1); off < A.length; off++) {
      var run = 0, start = -1;
      for (var i = 0; i <= B.length; i++) {
        var ai = off + i;
        var paired = (i < B.length && ai >= 0 && ai < A.length && A.charAt(ai) === B.charAt(i));
        if (paired) {
          if (run === 0) start = ai;
          run++;
        } else {
          if (run >= 3) {
            // B 已反向，故 A 与 B 上一致的片段即为互补配对段
            var seg = A.substr(start, run);
            var dG = run >= 2 ? nnDeltaG37(seg, opts) : 1.5;
            if (!isNaN(dG) && dG < best.dG) best = { dG: dG, run: run, offset: off };
          }
          run = 0;
        }
      }
    }
    return best;
  }

  /** 3′ 端是否参与二聚体（引物二聚体最危险的形态） */
  function dimerTouches3Prime(a, b, res) {
    if (!res || !res.run) return false;
    var A = String(a || '').toUpperCase();
    return (res.offset + res.run) >= A.length;   // 配对段延伸到 a 的 3′ 末端
  }

  /* =========================================================
   * 5. 分子量与浓度
   * ========================================================= */
  /**
   * 核酸分子量
   * @param {string} seq 序列
   * @param {string} kind 'ssDNA' | 'ssRNA' | 'dsDNA'
   * @param {object} [counts] 可选的碱基计数（baseCount 的结果）。
   *        传入可让调用方"一次计数、多次换算"，避免对同一条序列反复做 O(n) 统计。
   *        注意：counts 必须与 seq 的大小写状态一致（本函数内部会 toUpperCase）。
   */
  function mwNucleic(seq, kind, counts) {
    var s = seq.toUpperCase();
    var c = counts || baseCount(s);
    if (kind === 'ssRNA') {
      return c.A * RNA_RESIDUE.A + c.U * RNA_RESIDUE.U +
             c.C * RNA_RESIDUE.C + c.G * RNA_RESIDUE.G + WATER;
    }
    var t = c.T + c.U;
    var single = c.A * DNA_RESIDUE.A + t * DNA_RESIDUE.T +
                 c.C * DNA_RESIDUE.C + c.G * DNA_RESIDUE.G + WATER;
    if (kind === 'dsDNA') {
      // 输入为其中一条链：本链的 A+T 即双链中的 A:T 对数，G+C 即 G:C 对数
      // 单个 A:T 对 = 617.41 Da，G:C 对 = 618.39 Da（两条链残基质量之和）
      var nAT = c.A + t;
      var nGC = c.G + c.C;
      return nAT * 617.41 + nGC * 618.39;
    }
    return single;
  }

  function mwProtein(seq) {
    var total = WATER;
    for (var i = 0; i < seq.length; i++) {
      var r = AA_RESIDUE[String(seq).charAt(i).toUpperCase()];
      if (r !== undefined) total += r;
    }
    return total;
  }

  var OD_FACTOR = { dsDNA: 50, ssDNA: 33, ssRNA: 40 }; // µg/mL per 1.0 OD260

  /** 由 OD260 与稀释倍数换算浓度 */
  function fromOD(od, dilution, kind) {
    var f = OD_FACTOR[kind] || 50;
    return od * f * (dilution || 1); // µg/mL
  }

  /**
   * ng/µL → nM
   * 推导：ng/µL = µg/mL = 1e-3 g/L；mol/L = (ng/µL × 1e-3)/MW；nM = mol/L × 1e9
   * 即 nM = ng/µL × 1e6 / MW
   */
  function concNgUlToNm(ngPerUl, mw) {
    if (!mw) return 0;
    var molPerL = (ngPerUl / mw) * 1e-3; // (ng/µL)/(g/mol) → 10^-6 g/L 级别
    return molPerL * 1e9;
  }
  /** nM → ng/µL */
  function concNmToNgUl(nM, mw) {
    if (!mw) return 0;
    return (nM * 1e-9 * mw) * 1e6 * 1e-3;
  }

  /* =========================================================
   * 6. 限制性内切酶
   * ========================================================= */
  function iupacToRegex(site) {
    var re = '';
    for (var i = 0; i < site.length; i++) {
      var ch = site.charAt(i);
      re += IUPAC[ch] ? IUPAC[ch] : ch;
    }
    return new RegExp(re, 'g');
  }

  /**
   * 在序列中查找某酶的全部切点（支持重叠位点）
   * @param {boolean} circular 环状序列会额外识别跨越首尾接缝的位点
   *                  （把开头 site.length-1 个字符接到末尾再搜一次）
   */
  function findSites(seq, enzyme, circular) {
    var seqU = String(seq).toUpperCase();
    var site = enzyme.site.replace('^', '').toUpperCase();
    var cutOffset = enzyme.site.indexOf('^');
    if (cutOffset === -1) cutOffset = Math.floor(site.length / 2);

    var n = seqU.length;
    if (!n || !site.length) return [];
    // 环状时多接一小段开头序列，用于捕捉跨接缝位点；长度不超过 n-1 避免整序列重复匹配
    var tailPad = circular ? seqU.substr(0, Math.min(site.length - 1, Math.max(0, n - 1))) : '';
    var searchIn = seqU + tailPad;

    var re = iupacToRegex(site);
    var hits = [];
    var m;
    while ((m = re.exec(searchIn)) !== null) {
      if (m.index >= n && !circular) break;      // 线性序列不取跨接缝结果
      var wrap = m.index + site.length > n;      // 识别序列跨越了首尾接缝
      hits.push({
        start: (m.index % n) + 1,                        // 1-based 识别序列起点
        end: ((m.index + site.length - 1) % n) + 1,
        cut: (m.index + cutOffset) % n,                  // 切割发生在 cut 与 cut+1 之间 (0-based 边界)
        matched: m[0],
        wrap: wrap
      });
      re.lastIndex = m.index + 1;    // 允许重叠匹配
    }
    return hits;
  }

  /** 由切点集合计算片段长度 */
  function fragmentsFromCuts(len, cuts, circular) {
    var sorted = cuts.filter(function (p) { return p > 0 && p < len; })
                     .sort(function (a, b) { return a - b; });
    var uniq = [];
    for (var i = 0; i < sorted.length; i++) {
      if (i === 0 || sorted[i] !== sorted[i - 1]) uniq.push(sorted[i]);
    }
    if (uniq.length === 0) return [len];
    if (circular) {
      if (uniq.length === 1) return [len]; // 环上单切 → 线性化，一个全长片段
      var fr = [];
      for (var j = 1; j < uniq.length; j++) fr.push(uniq[j] - uniq[j - 1]);
      fr.push(len + uniq[0] - uniq[uniq.length - 1]);
      return fr.sort(function (a, b) { return b - a; });
    }
    var f = [uniq[0]];
    for (var k = 1; k < uniq.length; k++) f.push(uniq[k] - uniq[k - 1]);
    f.push(len - uniq[uniq.length - 1]);
    return f.sort(function (a, b) { return b - a; });
  }

  /**
   * 酶切反应
   * @param {string} seq
   * @param {Array} enzymes 酶对象数组 [{name, site}]
   * @param {boolean} circular 是否环状（质粒）
   */
  function digest(seq, enzymes, circular) {
    var len = seq.length;
    var allCuts = [];
    var perEnzyme = [];
    for (var i = 0; i < enzymes.length; i++) {
      var hits = findSites(seq, enzymes[i], circular);
      perEnzyme.push({
        enzyme: enzymes[i],
        hits: hits,
        count: hits.length
      });
      for (var j = 0; j < hits.length; j++) allCuts.push(hits[j].cut);
    }
    var frags = fragmentsFromCuts(len, allCuts, circular);
    return {
      length: len,
      circular: !!circular,
      perEnzyme: perEnzyme,
      fragments: frags
    };
  }

  /* =========================================================
   * 7. 引物评估
   * ========================================================= */
  function gcClamp(seq) {
    var s = seq.toUpperCase();
    var tail = s.slice(-5);
    var gc = (tail.match(/[GC]/g) || []).length;
    return { tail: tail, gc: gc };
  }

  /**
   * PCR 产物预测（多产物版本）。
   * 模板上引物可能有多个结合位点，会形成多种产物——诊断 PCR 里非特异条带
   * 往往就是短产物优先扩增的结果。因此这里返回全部组合，按产物长度升序
   * （越短越容易被优先扩增，越值得警惕），而不是只报最靠前的那个。
   * 环状模板上，反向引物位点位于正向引物"上游"同样能扩增出产物（绕接缝一圈），
   * 这类产物在质粒 PCR 中很常见，线性逻辑会直接漏掉。
   * @param {object} o {circular: 是否环状模板, maxProduct, maxResults}
   * @returns {Array} [{start,end,size,wrap,fwdPos,revPos,fwdIndex,revIndex}]
   */
  function pcrProducts(template, fwd, rev, o) {
    o = o || {};
    var circular = !!o.circular;
    var t = String(template || '').toUpperCase().replace(/[^ACGT]/g, '');
    var f = String(fwd || '').toUpperCase().replace(/[^ACGT]/g, '');
    var r = String(rev || '').toUpperCase().replace(/[^ACGT]/g, '');
    if (!t || !f || !r) return [];
    var maxProduct = o.maxProduct || 20000;
    var maxResults = o.maxResults || 60;
    var n = t.length;

    function findAll(hay, needle) {
      var res = [], i = hay.indexOf(needle);
      while (i !== -1) { res.push(i); i = hay.indexOf(needle, i + 1); }
      return res;
    }

    var rSite = revComp(r);
    var fSites = findAll(t, f);
    var rSites = findAll(t, rSite);
    var out = [];
    for (var i = 0; i < fSites.length; i++) {
      for (var j = 0; j < rSites.length; j++) {
        var fs = fSites[i], rs = rSites[j];
        var size, wrap;
        if (rs > fs) {
          // 常规：反向引物在正向引物下游
          size = rs + rSite.length - fs;
          wrap = false;
        } else if (circular && rs !== fs) {
          // 环状：产物从正向引物出发、越过序列末端绕回反向引物
          size = (n - fs) + rs + rSite.length;
          if (size > n) continue;        // 绕超过一整圈不成立
          wrap = true;
        } else {
          continue;                      // 线性模板不存在该产物
        }
        if (size > maxProduct) continue;
        out.push({
          start: fs + 1,
          end: rs + rSite.length,
          size: size,
          wrap: wrap,
          fwdPos: fs + 1,
          revPos: rs + 1,
          fwdIndex: i + 1,
          revIndex: j + 1
        });
      }
    }
    // 短产物在竞争中占优，排在前面
    out.sort(function (a, b) { return a.size - b.size; });
    return out.slice(0, maxResults);
  }

  /** 按产物记录取出产物序列（自动处理跨接缝） */
  function pcrProductSeq(template, p) {
    var t = String(template || '').toUpperCase().replace(/[^ACGT]/g, '');
    if (!p) return '';
    return p.wrap ? (t.slice(p.start - 1) + t.slice(0, p.end)) : t.slice(p.start - 1, p.end);
  }

  /* =========================================================
   * 8. 自动引物设计
   * ========================================================= */

  /**
   * 在模板的目标区域两端搜索候选引物并配对打分。
   * 过滤条件：长度、Tm（最近邻）、GC%、模板内唯一性、3′ 端 GC 夹
   * 排序依据：Tm 接近 60 °C、GC 接近 50%、发夹与二聚体越少越好、交叉二聚体越低越好
   * @param {string} template 模板序列
   * @param {object} o {start,end,minLen,maxLen,minTm,maxTm,minGc,maxGc,maxProduct,topN,searchWindow,conc,na,mg,dntp}
   * @returns {Array} 候选引物对，按 score 降序
   */
  function designPrimers(template, o) {
    o = o || {};
    var t = String(template || '').toUpperCase().replace(/U/g, 'T');
    var n = t.length;
    var minLen = o.minLen || 18, maxLen = o.maxLen || 27;
    var minTm = o.minTm || 55, maxTm = o.maxTm || 65;
    var minGc = o.minGc || 40, maxGc = o.maxGc || 60;
    var maxProduct = o.maxProduct || 4000;
    var topN = o.topN || 10;
    var searchWin = o.searchWindow || 80;
    var start = Math.max(0, o.start || 0);
    var end = Math.min(n, o.end === undefined ? n : o.end);
    if (end - start < 100 || !isValidNucleic(t)) return [];

    // 统计子串在模板中的出现次数（最多数到 3，用于判断唯一性）
    function occurrences(sub) {
      var c = 0, i = t.indexOf(sub);
      while (i !== -1 && c < 3) { c++; i = t.indexOf(sub, i + 1); }
      return c;
    }

    /**
     * 结合位点周边的局部 GC（前后各 50 bp）。
     * 模板局部 GC 过高（二级结构、难解链）或过低（易非特异结合）都会降低扩增效率，
     * 这里作为"上下文"指标参与打分，而不只是看引物自身的 GC%。
     */
    function localGc(pos, len) {
      var a = Math.max(0, pos - 50);
      var b = Math.min(n, pos + len + 50);
      if (b <= a) return 50;
      return gcContent(t.substr(a, b - a));
    }

    function evaluate(seq, isRev, pos) {
      if (!/^[ACGT]+$/.test(seq)) return null;
      var gc = gcContent(seq);
      if (gc < minGc || gc > maxGc) return null;
      var tm = tmNearestNeighbor(seq, o);
      if (isNaN(tm) || tm < minTm || tm > maxTm) return null;

      var clamp = gcClamp(seq);
      var hp = hairpinDeltaG(seq, o);
      var sd = dimerDeltaG(seq, seq, o);
      var occ = occurrences(isRev ? revComp(seq) : seq);
      var lgc = localGc(pos, seq.length);

      // ΔG 为负值（越负越不稳定），直接以负系数折进评分
      var score = 100
        - Math.abs(tm - 60) * 2
        - Math.abs(gc - 50) * 0.8
        + hp.dG * 2.5
        + sd.dG * 2.0;
      if (clamp.gc < 2) score -= 8;
      else if (clamp.gc >= 3) score += 4;
      var last = seq.charAt(seq.length - 1);
      if (last === 'G' || last === 'C') score += 3;
      if (occ > 1) score -= 30;                 // 模板中存在多个结合位点 → 非特异风险
      if (lgc > 70) score -= (lgc - 70) * 0.6;  // 局部 GC 过高，二级结构风险
      else if (lgc < 30) score -= (30 - lgc) * 0.6;

      return {
        seq: seq, gc: gc, tm: tm, clamp: clamp.gc,
        hairpinDg: hp.dG, selfDimerDg: sd.dG, occ: occ,
        localGc: lgc, score: score
      };
    }

    var fwd = [], rev = [], i, j, e;
    // 正向：目标区域起点附近，沿正链取
    for (i = start; i <= start + searchWin && i + minLen <= n; i++) {
      for (j = minLen; j <= maxLen && i + j <= n; j++) {
        e = evaluate(t.substr(i, j), false, i);
        if (e) { e.pos = i; fwd.push(e); }
      }
    }
    // 反向：目标区域终点附近，取互补链的反向互补
    for (i = end; i >= end - searchWin && i - minLen >= 0; i--) {
      for (j = minLen; j <= maxLen && i - j >= 0; j++) {
        e = evaluate(revComp(t.substr(i - j, j)), true, i - j);
        if (e) { e.pos = i; rev.push(e); }
      }
    }
    if (!fwd.length || !rev.length) return [];

    fwd.sort(function (a, b) { return b.score - a.score; });
    rev.sort(function (a, b) { return b.score - a.score; });
    fwd = fwd.slice(0, 40);
    rev = rev.slice(0, 40);

    var pairs = [];
    for (var a = 0; a < fwd.length; a++) {
      for (var b = 0; b < rev.length; b++) {
        var f = fwd[a], r = rev[b];
        var size = r.pos - f.pos;
        if (size < 60 || size > maxProduct) continue;
        var tmDiff = Math.abs(f.tm - r.tm);
        if (tmDiff > 3) continue;
        var cross = dimerDeltaG(f.seq, r.seq, o);
        var crossR = dimerDeltaG(r.seq, f.seq, o);
        var t3f = dimerTouches3Prime(f.seq, r.seq, cross);
        var t3r = dimerTouches3Prime(r.seq, f.seq, crossR);
        var worst3 = t3f ? cross.dG : (t3r ? crossR.dG : 0);
        pairs.push({
          fwd: f, rev: r, size: size,
          start: f.pos + 1, end: r.pos,
          tmDiff: tmDiff,
          crossDg: cross.dG, crossRun: cross.run,
          touch3: !!(t3f || t3r), dg3: worst3,
          score: f.score + r.score + cross.dG * 2.5 + worst3 * 3
        });
      }
    }
    pairs.sort(function (x, y) { return y.score - x.score; });
    return pairs.slice(0, topN);
  }

  /* =========================================================
   * 9. 双序列全局比对（Needleman-Wunsch）
   * ========================================================= */

  /**
   * Needleman-Wunsch 全局比对（线性 gap 罚分）。
   * 复杂度 O(n·m) 时间与空间；超过 maxCells 个矩阵单元时拒绝执行，避免卡死浏览器。
   * 完整多序列比对（MSA）需要渐进式比对或迭代策略，规模与复杂度都是另一个量级，
   * 不在本工具范围内。
   * @returns {{score,identity,gaps,alignedA,alignedB,midline,length,match,lenA,lenB}|null|{tooBig:true}}
   */
  function alignPair(a, b, opts) {
    opts = opts || {};
    var A = String(a || '').toUpperCase().replace(/[^A-Z]/g, '');
    var B = String(b || '').toUpperCase().replace(/[^A-Z]/g, '');
    var n = A.length, m = B.length;
    if (!n || !m) return null;

    var MATCH = opts.match === undefined ? 2 : opts.match;
    var MIS = opts.mismatch === undefined ? -1 : opts.mismatch;
    var GAP = opts.gap === undefined ? -2 : opts.gap;
    var maxCells = opts.maxCells || 4000000;   // 约 16 MB（Int32），够 2000×2000
    if (n * m > maxCells) return { tooBig: true, n: n, m: m, cells: n * m, maxCells: maxCells };

    var cols = m + 1;
    var F = new Int32Array((n + 1) * cols);
    var TB = new Uint8Array((n + 1) * cols);   // 0=对角 1=上 2=左
    var i, j;

    for (j = 1; j <= m; j++) { F[j] = j * GAP; TB[j] = 2; }
    for (i = 1; i <= n; i++) { F[i * cols] = i * GAP; TB[i * cols] = 1; }

    for (i = 1; i <= n; i++) {
      var ai = A.charAt(i - 1);
      for (j = 1; j <= m; j++) {
        var diag = F[(i - 1) * cols + (j - 1)] + (ai === B.charAt(j - 1) ? MATCH : MIS);
        var up = F[(i - 1) * cols + j] + GAP;
        var left = F[i * cols + (j - 1)] + GAP;
        var best = diag, tb = 0;
        if (up > best) { best = up; tb = 1; }
        if (left > best) { best = left; tb = 2; }
        F[i * cols + j] = best;
        TB[i * cols + j] = tb;
      }
    }

    // 回溯
    var oa = [], ob = [], om = [];
    i = n; j = m;
    var ident = 0, gaps = 0, alnLen = 0;
    while (i > 0 || j > 0) {
      var t = (i === 0) ? 2 : (j === 0 ? 1 : TB[i * cols + j]);
      if (i > 0 && j > 0 && t === 0) {
        var ca = A.charAt(i - 1), cb = B.charAt(j - 1);
        oa.push(ca); ob.push(cb);
        if (ca === cb) { om.push('|'); ident++; } else om.push('.');
        i--; j--;
      } else if (t === 1 && i > 0) {
        oa.push(A.charAt(i - 1)); ob.push('-'); om.push(' '); gaps++; i--;
      } else {
        oa.push('-'); ob.push(B.charAt(j - 1)); om.push(' '); gaps++; j--;
      }
      alnLen++;
    }
    oa.reverse(); ob.reverse(); om.reverse();

    return {
      score: F[n * cols + m],
      identity: alnLen ? (ident / alnLen) * 100 : 0,
      gaps: gaps,
      match: ident,
      length: alnLen,
      lenA: n, lenB: m,
      alignedA: oa.join(''),
      alignedB: ob.join(''),
      midline: om.join('')
    };
  }

  /* =========================================================
   * 10. 虚拟克隆（酶切 → 片段选择 → 连接）
   * ========================================================= */

  /**
   * 某个酶切点产生的末端形态。
   * 从识别序列里 `^` 的位置就能推出：切点左右两侧长度不等 → 产生突出端，
   * 下游更长是 5′ 突出（如 EcoRI G^AATTC → 5′-AATT），上游更长是 3′ 突出
   * （如 PstI CTGCA^G → 3′-TGCA），等长则平端（如 SmaI CCC^GGG）。
   */
  function enzymeEnds(enzyme) {
    var site = String(enzyme.site).replace('^', '').toUpperCase();
    var off = String(enzyme.site).indexOf('^');
    if (off === -1) off = Math.floor(site.length / 2);
    var upstream = site.slice(0, off);
    var downstream = site.slice(off);
    var olen = Math.abs(upstream.length - downstream.length);
    var type, overhang = '';
    if (!olen) {
      type = 'blunt';
    } else if (upstream.length < downstream.length) {
      type = '5';
      overhang = downstream.slice(0, olen);
    } else {
      type = '3';
      overhang = upstream.slice(upstream.length - olen);
    }
    return {
      type: type,
      length: olen,
      overhang: overhang,
      // 简并位点（如 AvaI C^YCGRG）推不出确切的突出序列，连接判断需谨慎
      degenerate: /[^ACGT]/.test(overhang),
      enzyme: enzyme.name || ''
    };
  }

  /**
   * 两个末端能否连接：
   *  - 双方都是平端 → 可以（任何平端之间都能连）
   *  - 同型（都是 5′ 或都是 3′）且突出序列相同 → 可以
   *    （同一个酶切出的两端天然自洽；BamHI/BglII 这类同尾酶也属于此列）
   *  - 类型不同或突出序列不同 → 不能直接连
   */
  function endsCompatible(a, b) {
    if (!a || !b) return false;
    if (a.degenerate || b.degenerate) return false;   // 简并突出端无法判定
    if (a.type === 'blunt' && b.type === 'blunt') return true;
    if (a.type !== b.type) return false;
    return a.overhang === b.overhang;
  }

  /**
   * 单酶切并返回带末端信息的片段。
   * 环状时首尾相接成环，因此第 1 个片段的左端来自最后一个切点。
   * @returns {Array} [{ seq, length, leftCut, rightCut, leftEnd, rightEnd, wrap }]
   */
  function digestOneDetailed(seq, enzyme, circular) {
    var n = seq.length;
    if (!n) return [];
    var hits = findSites(seq, enzyme, circular);
    var cuts = [];
    hits.forEach(function (h) { if (h.cut > 0 && h.cut < n) cuts.push(h.cut); });
    cuts = cuts.filter(function (v, i, a) { return a.indexOf(v) === i; }).sort(function (a, b) { return a - b; });
    var ends = enzymeEnds(enzyme);

    if (!cuts.length) {
      // 未切开：环状是完整的环，线性是一条完整序列
      return [{ seq: seq, length: n, leftCut: 0, rightCut: n, leftEnd: null, rightEnd: null, uncut: true }];
    }
    var out = [];
    if (!circular) {
      var prev = 0;
      cuts.forEach(function (c) {
        out.push({
          seq: seq.slice(prev, c), length: c - prev,
          leftCut: prev, rightCut: c,
          leftEnd: prev === 0 ? null : ends, rightEnd: ends
        });
        prev = c;
      });
      out.push({
        seq: seq.slice(prev), length: n - prev,
        leftCut: prev, rightCut: n,
        leftEnd: ends, rightEnd: null
      });
    } else {
      for (var i = 0; i < cuts.length; i++) {
        var a = cuts[i], b = cuts[(i + 1) % cuts.length];
        var sub, wrap = false;
        if (b > a) sub = seq.slice(a, b);
        else { sub = seq.slice(a) + seq.slice(0, b); wrap = true; }   // 跨接缝的那一段
        out.push({
          seq: sub, length: sub.length,
          leftCut: a, rightCut: b,
          leftEnd: ends, rightEnd: ends, wrap: wrap
        });
      }
    }
    return out;
  }

  /**
   * 虚拟克隆：把载体骨架上的一段替换为插入片段（最常见的定点克隆）。
   * 骨架用两个酶切除掉中间一段，插入片段两端也由酶切产生，
   * 这里检查两个连接点的黏端/平端是否匹配，并给出重组后的序列。
   * @param {object} backbone { seq, leftEnd, rightEnd } 骨架自己首尾的连接端
   * @param {object} insert { seq, leftEnd, rightEnd }
   * @returns {{ok, seq, size, junctions, reason}}
   */
  function ligateInsert(backbone, insert) {
    if (!backbone || !insert || !backbone.seq || !insert.seq) {
      return { ok: false, reason: '缺少载体骨架或插入片段的序列' };
    }
    var junctions = [];

    // 连接点 1：骨架左端 ← 插入右端
    var j1 = endsCompatible(insert.rightEnd, backbone.leftEnd);
    junctions.push({
      name: '插入 3′ ↔ 骨架 5′',
      byInsert: insert.rightEnd, byBackbone: backbone.leftEnd, ok: j1
    });

    // 连接点 2：插入左端 ← 骨架右端
    var j2 = endsCompatible(backbone.rightEnd, insert.leftEnd);
    junctions.push({
      name: '骨架 3′ ↔ 插入 5′',
      byBackbone: backbone.rightEnd, byInsert: insert.leftEnd, ok: j2
    });

    var allOk = junctions.every(function (j) { return j.ok; });
    var seq = backbone.seq + insert.seq;
    return {
      ok: allOk,
      seq: seq,
      size: seq.length,
      junctions: junctions,
      reason: allOk ? '' : '连接点末端不匹配，无法直接连接（可考虑补平、加接头或换用同尾酶）'
    };
  }

  /**
   * 多序列比对（中心星 / center-star 算法）。
   *
   * 思路：
   *  1. 两两比对算距离，选与其他序列平均距离最小的那条作「中心」；
   *  2. 中心依次与其余序列两两比对；
   *  3. 把每条比对投影到统一坐标系：中心每个位置之前需要多少 gap 列，
   *     取所有序列的逐位最大值，再把各序列摆进去。
   *
   * gap 列的处置规则：由中心引入的 gap 列，所有序列在该列都是 `-`；
   * 中心有字符的位置，某条序列若在该处缺碱基则记 `-`。
   * 这样能保证所有行等长、去 gap 后都能还原成原始序列。
   *
   * @returns {{rows, count, centers, length, pairwise}|{tooBig:true}}
   */
  function alignMultiple(seqs, opts) {
    opts = opts || {};
    var list = (seqs || []).map(function (s) {
      return String(s || '').toUpperCase().replace(/[^A-Z]/g, '');
    }).filter(function (s) { return s.length > 0; });
    if (list.length < 2) return null;
    var maxSeqs = opts.maxSeqs || 20;
    if (list.length > maxSeqs) return { tooBig: true, reason: '序列条数超过 ' + maxSeqs };

    var i, j, n = list.length;

    // 两两比对：既用于选中心，也用于最后的相似度展示
    var pairIds = [];
    for (i = 0; i < n; i++) {
      for (j = i + 1; j < n; j++) {
        var p = alignPair(list[i], list[j], opts);
        if (!p) continue;
        if (p.tooBig) return { tooBig: true, cells: p.cells, maxCells: p.maxCells };
        pairIds.push({ a: i, b: j, identity: p.identity, gaps: p.gaps });
      }
    }

    // 选中心：与其余序列相似度之和最大（即距离和最小）
    var center = 0, best = -Infinity;
    for (i = 0; i < n; i++) {
      var sum = 0, cnt = 0;
      pairIds.forEach(function (pr) {
        if (pr.a === i || pr.b === i) { sum += pr.identity; cnt++; }
      });
      var avg = cnt ? sum / cnt : 0;
      if (avg > best) { best = avg; center = i; }
    }
    var centerSeq = list[center];

    // 每条序列在中心坐标系下的展开：gapBefore[p] = 位置 p 之前的 gap 数，
    // chars[p] = 位置 p 对应的字符（可能本身就是 '-'）
    var expansions = [];
    for (i = 0; i < n; i++) {
      if (i === center) continue;
      var al = alignPair(centerSeq, list[i], opts);
      if (!al || al.tooBig) continue;
      // gapBefore[p] = 中心位置 p 之前的那几列，逐列记录本条序列的字符。
      // 关键：这些列里本条序列可能是插入的碱基（不是 '-'），
      // 早期实现只数了 gap 个数、丢掉了插入碱基，导致去 gap 后还原不出原序列。
      var gapBefore = [], chars = [], cur = [];
      for (var q = 0; q < al.alignedA.length; q++) {
        if (al.alignedA.charAt(q) === '-') {
          cur.push(al.alignedB.charAt(q));
        } else {
          gapBefore.push(cur); cur = [];
          chars.push(al.alignedB.charAt(q));
        }
      }
      gapBefore.push(cur);
      expansions.push({ idx: i, gapBefore: gapBefore, chars: chars });
    }

    // 逐位取 gap 列数的最大值
    var maxGap = new Array(centerSeq.length + 1);
    for (i = 0; i <= centerSeq.length; i++) {
      var mx = 0;
      expansions.forEach(function (e) {
        var len = (e.gapBefore[i] || []).length;
        if (len > mx) mx = len;
      });
      maxGap[i] = mx;
    }

    function buildRow(exp) {
      var out = '';
      for (var p = 0; p < centerSeq.length; p++) {
        var cols = exp ? (exp.gapBefore[p] || []) : [];
        for (var g = 0; g < maxGap[p]; g++) {
          out += (g < cols.length) ? cols[g] : '-';
        }
        out += exp ? exp.chars[p] : centerSeq.charAt(p);
      }
      // 末尾的 gap 列
      var tail = exp ? (exp.gapBefore[centerSeq.length] || []) : [];
      for (var g2 = 0; g2 < maxGap[centerSeq.length]; g2++) {
        out += (g2 < tail.length) ? tail[g2] : '-';
      }
      return out;
    }

    var rows = new Array(n);
    rows[center] = buildRow(null);
    expansions.forEach(function (e) { rows[e.idx] = buildRow(e); });

    // 兜底对齐：拼接过程中若出现长度差异，右侧补 gap
    var width = 0;
    rows.forEach(function (r) { if (r && r.length > width) width = r.length; });
    rows = rows.map(function (r) {
      if (r === undefined || r === null) return new Array(width + 1).join('-');
      while (r.length < width) r += '-';
      return r;
    });

    // 每个位点的一致度（用于着色）
    var matchCols = 0;
    for (var c = 0; c < width; c++) {
      var first = rows[0].charAt(c), same = true;
      for (var k = 1; k < rows.length; k++) if (rows[k].charAt(c) !== first) { same = false; break; }
      if (same && first !== '-') matchCols++;
    }

    return {
      rows: rows,
      count: n,
      center: center,
      centerIndexInRows: center,
      length: width,
      columns: width,
      conserved: matchCols,
      consensusPct: width ? (matchCols / width) * 100 : 0,
      pairwise: pairIds
    };
  }

  /* =========================================================
   * 工具函数
   * ========================================================= */
  /* =========================================================
   * 10. 蛋白理化性质（pI / GRAVY / 不稳定指数 / 脂肪族指数 / 消光系数）
   *     参数出处：EMBOSS pKa；Kyte-Doolittle 1982；Ikai 1980；
   *     Guruprasad 1990 DIWV；Pace 1995 消光系数
   * ========================================================= */
  var AA_ORDER = 'ACDEFGHIKLMNPQRSTVWY';
  var KD_HYDRO = {
    A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5, Q: -3.5, E: -3.5, G: -0.4,
    H: -3.2, I: 4.5, L: 3.8, K: -3.9, M: 1.9, F: 2.8, P: -1.6, S: -0.8,
    T: -0.7, W: -0.9, Y: -1.3, V: 4.2
  };
  // EMBOSS 采用的 Bjellqvist pKa（侧链 + 末端）
  var PKA = { NTER: 8.6, CTER: 3.6, K: 10.8, R: 12.5, H: 6.5, D: 3.9, E: 4.1, C: 8.5, Y: 10.1 };
  // Guruprasad et al. 1990 二肽不稳定权重（行=第一个残基，列顺序 ACDEFGHIKLMNPQRSTVWY）
  var DIWV_ROWS = {
    A: [1, 44.94, -7.49, 1, 1, 1, -7.49, 1, 1, 1, 1, 1, 20.26, 1, 1, 1, 1, 1, 1, 1],
    C: [1, 1, 20.26, 1, 1, 1, 33.60, 1, 1, 20.26, 33.60, 1, 20.26, -6.54, 1, 1, 33.60, -6.54, 24.68, 1],
    D: [1, 1, 1, 1, -6.54, 1, 1, 1, -7.49, 1, 1, 1, 1, 1, -6.54, 20.26, -14.03, 1, 1, 1],
    E: [1, 44.94, 20.26, 33.60, 1, 1, -6.54, 20.26, 1, 1, 1, 1, 20.26, 20.26, 1, 20.26, 1, 1, -14.03, 1],
    F: [1, 1, 13.34, 1, 1, 1, 1, 1, -14.03, 1, 1, 1, 20.26, 1, 1, 1, 1, 1, 1, 33.601],
    G: [-7.49, 1, 1, -6.54, 1, 13.34, 1, -7.49, -7.49, 1, 1, -7.49, 1, 1, 1, 1, -7.49, 1, 13.34, -7.49],
    H: [1, 1, 1, 1, -9.37, -9.37, 1, 44.94, 24.68, 1, 1, 24.68, -1.88, 1, 1, 1, -6.54, 1, -1.88, 44.94],
    I: [1, 1, 1, 44.94, 1, 1, 13.34, 1, -7.49, 20.26, 1, 1, -1.88, 1, 1, 1, 1, -7.49, 1, 1],
    K: [1, 1, 1, 1, 1, -7.49, 1, 1, 1, -7.49, 33.60, 1, -6.54, 24.64, 33.60, 1, 1, -7.49, 1, 1],
    L: [1, 1, 1, 1, 1, 1, 1, 1, -7.49, 1, 1, 1, 20.26, 33.60, 1, 1, 24.68, 1, 1, 1],
    M: [13.34, 1, 1, 1, 1, 1, 58.28, 1, 1, 1, -1.88, 1, 44.94, -6.54, -6.54, 44.94, -1.88, 1, 1, 24.68],
    N: [1, -1.88, 1, 1, -14.03, -14.03, 1, 44.94, 24.68, 1, 1, 1, -1.88, -6.54, 1, 1, 1, 1, -9.37, 1],
    P: [20.26, -6.54, -6.54, 18.38, 20.26, 1, 1, 1, 1, 1, -6.54, 1, 20.26, 1, -6.54, 20.26, 1, 20.26, -1.88, 1],
    Q: [1, -6.54, 20.26, 20.26, -6.54, 1, 1, 1, 1, 1, 1, 1, 20.26, 20.26, 1, 44.94, 1, 1, -6.54, -6.54],
    R: [1, 1, 1, 1, 1, -7.49, 20.26, 1, 1, 1, 13.34, 13.34, 20.26, 20.26, 58.28, 44.94, 1, 1, 58.28, -6.54],
    S: [1, 33.60, 1, 20.26, 1, 1, 1, 1, 1, 1, 1, 1, 44.94, 20.26, 1, 20.26, 1, 1, 1, 1],
    T: [1, 1, 1, 20.26, 13.34, -7.49, 1, 1, 1, 1, 1, -14.03, 1, -6.54, 1, 1, 1, 1, -14.03, 1],
    V: [1, 1, -14.03, 1, 1, -7.49, 1, 1, -1.88, 1, 1, 1, 20.26, 1, 1, 1, 1, 1, 1, -6.54],
    W: [-14.03, 1, 1, 1, 1, -9.37, 24.68, 1, 24.68, 13.34, 24.68, 13.34, 1, 1, 1, 1, -14.03, 1, 1, -7.49],
    Y: [24.68, 1, 24.68, -6.54, 1, -7.49, 13.34, 1, 44.94, 1, 1, 1, 13.34, 1, -15.91, 1, -7.49, 1, -9.37, 13.34]
  };
  var DIWV = {};
  (function () {
    AA_ORDER.split('').forEach(function (a, i) {
      DIWV[a] = {};
      AA_ORDER.split('').forEach(function (b, j) { DIWV[a][b] = DIWV_ROWS[a][j]; });
    });
  })();

  function aaCounts(aa) {
    var c = {};
    String(aa).toUpperCase().split('').forEach(function (r) {
      if (AA_ORDER.indexOf(r) >= 0) c[r] = (c[r] || 0) + 1;
    });
    return c;
  }

  /** 指定 pH 下蛋白净电荷（Henderson-Hasselbalch） */
  function proteinNetCharge(counts, len, pH) {
    var pos = 1 / (1 + Math.pow(10, pH - PKA.NTER)); // N 端
    var neg = 1 / (1 + Math.pow(10, PKA.CTER - pH)); // C 端
    ['K', 'R', 'H'].forEach(function (r) {
      pos += (counts[r] || 0) / (1 + Math.pow(10, pH - PKA[r]));
    });
    ['D', 'E', 'C', 'Y'].forEach(function (r) {
      neg += (counts[r] || 0) / (1 + Math.pow(10, PKA[r] - pH));
    });
    return pos - neg;
  }

  /** 等电点：在 0~14 内二分搜索净电荷零点 */
  function proteinPI(counts, len) {
    if (len === 0) return null;
    var lo = 0, hi = 14;
    if (proteinNetCharge(counts, len, 14) > 0) return 14;
    if (proteinNetCharge(counts, len, 0) < 0) return 0;
    for (var i = 0; i < 80; i++) {
      var mid = (lo + hi) / 2;
      if (proteinNetCharge(counts, len, mid) > 0) lo = mid; else hi = mid;
    }
    return (lo + hi) / 2;
  }

  /**
   * 蛋白理化性质汇总。输入氨基酸序列（允许 * 与非法字符，自动忽略）。
   * 返回 null 表示没有任何有效残基。
   */
  function proteinProperties(aa) {
    // 只保留 20 种标准残基：*、数字、X/Z/B/J 等非常规字母一律剔除，避免后续查表越界
    var seq = String(aa).toUpperCase().split('').filter(function (r) {
      return AA_ORDER.indexOf(r) >= 0;
    }).join('');
    var len = seq.length;
    if (!len) return null;
    var counts = aaCounts(seq);
    var mw = mwProtein(seq);
    var pI = proteinPI(counts, len);
    // GRAVY：Kyte-Doolittle 疏水性平均
    var gravy = 0;
    seq.split('').forEach(function (r) { gravy += (KD_HYDRO[r] || 0); });
    gravy /= len;
    // 脂肪族指数（Ikai 1980）：X(A) + 2.9·X(V) + 3.9·[X(I)+X(L)]，X 为摩尔百分数
    var pct = function (r) { return 100 * (counts[r] || 0) / len; };
    var aliphatic = pct('A') + 2.9 * pct('V') + 3.9 * (pct('I') + pct('L'));
    // 不稳定指数（Guruprasad 1990）：10/L × Σ DIWV(xi, xi+1)
    var di = 0;
    for (var i = 0; i < len - 1; i++) di += DIWV[seq[i]][seq[i + 1]];
    var instability = len > 1 ? (10 / len) * di : 0;
    // 消光系数（Pace 1995）：ε280 = 5500·W + 1490·Y (+ 125·C 若全部形成二硫键)
    var extRed = 5500 * (counts.W || 0) + 1490 * (counts.Y || 0);
    var extOx = extRed + 125 * (counts.C || 0);
    return {
      length: len,
      mw: mw,
      pI: pI,
      charge7: proteinNetCharge(counts, len, 7),
      gravy: gravy,
      aliphatic: aliphatic,
      instability: instability,
      unstable: instability > 40,
      extinctionOx: extOx,
      extinctionRed: extRed,
      abs01Ox: extOx / mw,   // 1 g/L 溶液的 A280
      abs01Red: extRed / mw,
      counts: counts
    };
  }

  /* =========================================================
   * 11. 密码子优化 / 反向翻译（E. coli K12 高表达偏好）
   * ========================================================= */
  var ECOLI_PREF = {
    F: 'TTC', L: 'CTG', I: 'ATC', M: 'ATG', V: 'GTG', S: 'AGC', P: 'CCG',
    T: 'ACC', A: 'GCG', Y: 'TAC', H: 'CAC', Q: 'CAG', N: 'AAC', K: 'AAA',
    D: 'GAC', E: 'GAA', C: 'TGC', W: 'TGG', R: 'CGT', G: 'GGC', '*': 'TAA'
  };

  /**
   * 把蛋白序列反向翻译为 E. coli 偏好密码子的 DNA。
   * 返回 { dna, skipped, backTranslateOk }；跳过非法残基并在 skipped 中计数。
   */
  function reverseTranslate(aa) {
    var seq = String(aa).toUpperCase().replace(/[^A-Z*]/g, '');
    var dna = '', skipped = 0;
    for (var i = 0; i < seq.length; i++) {
      var codon = ECOLI_PREF[seq[i]];
      if (codon) dna += codon;
      else skipped++;
    }
    var back = translate(dna).replace(/\*$/, '');
    var origClean = seq.replace(/\*+$/, '');
    return { dna: dna, skipped: skipped, backTranslateOk: back === origClean };
  }

  /* =========================================================
   * 12. 连接反应摩尔比计算
   * ========================================================= */
  /**
   * 载体 / 插入片段用量。约定双链 DNA 平均 660 g/mol/bp（实验室通用系数）。
   * @param ratio 插入 : 载体 摩尔比（如 3）
   * @returns {vectorPmol, insertPmol, insertNg}
   */
  function ligationAmounts(vectorBp, insertBp, vectorNg, ratio) {
    vectorBp = Number(vectorBp); insertBp = Number(insertBp);
    vectorNg = Number(vectorNg); ratio = Number(ratio);
    if (!(vectorBp > 0) || !(insertBp > 0) || !(vectorNg > 0) || !(ratio > 0)) return null;
    var vectorPmol = vectorNg * 1000 / (vectorBp * 660);
    var insertPmol = vectorPmol * ratio;
    return {
      vectorPmol: vectorPmol,
      insertPmol: insertPmol,
      insertNg: insertPmol * (insertBp * 660) / 1000
    };
  }

  /* =========================================================
   * 13. 定点突变引物设计（Q5® SDM 型：背对背、突变在正向 5′ 端）
   * ========================================================= */
  /**
   * @param template 模板序列（质粒时勾选 circular）
   * @param pos      突变起点（1-based，指向被替换区第一个碱基）
   * @param newBases 替换为的碱基（'' = 缺失；长度≠oldLen 时为插入/替换）
   * @param opts     { oldLen=1, targetTm=62, minAnneal=18, maxAnneal=30, circular }
   */
  function designMutPrimers(template, pos, newBases, opts) {
    opts = opts || {};
    var oldLen = opts.oldLen === undefined ? 1 : Math.max(0, opts.oldLen | 0);
    var targetTm = opts.targetTm || 62;
    var minAnneal = opts.minAnneal || 18;
    var maxAnneal = opts.maxAnneal || 30;
    var seq = String(template).toUpperCase().replace(/[^ACGTU]/g, '').replace(/U/g, 'T');
    var n = seq.length;
    newBases = String(newBases || '').toUpperCase().replace(/[^ACGT]/g, '');
    if (!n) return { ok: false, reason: '模板为空' };
    pos = Math.round(Number(pos));
    if (!(pos >= 1 && pos <= n)) return { ok: false, reason: '突变位置超出模板范围（1–' + n + '）' };
    if (oldLen > n) return { ok: false, reason: '被替换碱基数（' + oldLen + '）超过模板全长（' + n + '）' };
    if (pos + oldLen - 1 > n && !opts.circular) return { ok: false, reason: '被替换区超出模板末端' };
    var doubled = seq + seq; // 环状取模索引
    function at(i) { return doubled.charAt(((i % n) + n) % n); }
    function sub(from0, len) { // 0-based 起点、可跨接缝
      var s = '';
      for (var i = 0; i < len; i++) s += at(from0 + i);
      return s;
    }
    // 正向：突变碱基置于 5′ 端，下游退火区延伸到 Tm 达标
    var tmOpts = { dnaConc: (opts.dnaConc || 500), naConc: (opts.naConc || 50) };
    function extendAnneal(build) {
      var best = null;
      for (var L = minAnneal; L <= maxAnneal; L++) {
        var anneal = build(L);
        var tm = tmNearestNeighbor(anneal, tmOpts);
        best = { anneal: anneal, tm: tm, len: L };
        if (tm >= targetTm) break;
      }
      return best;
    }
    var fwdA = extendAnneal(function (L) { return sub(pos - 1 + oldLen, L); });
    var revA = extendAnneal(function (L) { return revComp(sub(pos - 1 - L, L)); });
    if (!fwdA || !revA) return { ok: false, reason: '无法构造退火区' };
    var fwd = newBases + fwdA.anneal;
    var rev = revA.anneal;
    var warnings = [];
    if (fwdA.tm < targetTm) warnings.push('正向退火区已达 ' + maxAnneal + ' nt 上限，Tm 仍只有 ' + fwdA.tm.toFixed(1) + ' °C（AT 富集区？建议换 Q5 等长延伸酶并降低退火温度）');
    if (revA.tm < targetTm) warnings.push('反向退火区已达 ' + maxAnneal + ' nt 上限，Tm 只有 ' + revA.tm.toFixed(1) + ' °C');
    if (fwd.length > 60) warnings.push('正向引物总长 ' + fwd.length + ' nt 超过 60，合成成本与错误率上升，考虑把突变拆到反向引物上');
    var mutDesc;
    if (oldLen === 1 && newBases.length === 1) mutDesc = seq.charAt(pos - 1) + pos + newBases;
    else if (newBases.length === 0) mutDesc = '缺失 ' + pos + '–' + (pos + oldLen - 1) + '（' + oldLen + ' bp）';
    else if (oldLen === 0) mutDesc = '在 ' + (pos - 1) + ' 与 ' + pos + ' 之间插入 ' + newBases.length + ' bp';
    else mutDesc = '替换 ' + pos + '–' + (pos + oldLen - 1) + ' 为 ' + newBases.length + ' bp';
    return {
      ok: true,
      fwd: fwd, rev: rev,
      fwdTail: newBases, fwdAnneal: fwdA.anneal, revAnneal: revA.anneal,
      tmF: fwdA.tm, tmR: revA.tm,
      productLen: n - oldLen + newBases.length,
      mutDesc: mutDesc,
      warnings: warnings
    };
  }

  /* =========================================================
   * 14. Motif / IUPAC 模式搜索
   * ========================================================= */
  /**
   * 在序列中查找 IUPAC 简并 motif（含重叠匹配）。
   * @param opts { circular, bothStrands, maxResults=50000 }
   *        负链命中会映射回正链坐标；跨接缝命中带 wrap 标记。
   * @returns [{ pos(1-based), strand:'+'|'-', match, length, wrap }]
   *          数组上另有非枚举属性 `truncated`，为 true 表示命中数触达上限被截断。
   */
  function findMotifs(seq, pattern, opts) {
    opts = opts || {};
    var seqU = String(seq).toUpperCase();
    var pat = String(pattern).toUpperCase().replace(/\s/g, '');
    if (!seqU || !pat) return [];
    var maxResults = opts.maxResults === undefined ? 50000 : opts.maxResults;
    if (!(maxResults > 0)) maxResults = 50000;
    var re = iupacToRegex(pat);
    var m = pat.length;
    var n = seqU.length;
    var out = [];
    var truncated = false;
    var N;
    function scan(hay, strand) {
      N = hay.length;
      re.lastIndex = 0;
      var hit;
      while ((hit = re.exec(hay)) !== null) {
        // 正链：扩展区里的命中是 seq 已有命中的重复，遇到即止
        if (strand === '+' && opts.circular && hit.index >= n) break;
        var len = hit[0].length;
        var pos1;
        if (strand === '+') {
          pos1 = hit.index + 1;
        } else {
          // 负链命中覆盖 hay 的 [N-i-len, N-1-i]，起点映射回正链坐标；
          // 取模使其在跨接缝时同样成立（早期写法 n-i-len+1 在环状跨接缝会算出 0 而丢命中）
          pos1 = (((N - hit.index - len) % n) + n) % n + 1;
        }
        var wrap = !!(opts.circular && pos1 + len - 1 > n);
        if (pos1 >= 1 && pos1 <= n) {
          if (out.length >= maxResults) { truncated = true; return; }
          out.push({ pos: pos1, strand: strand, match: hit[0], length: len, wrap: wrap });
        }
        re.lastIndex = hit.index + 1; // 允许重叠
      }
    }
    var hay = opts.circular ? seqU + seqU.slice(0, m - 1) : seqU;
    scan(hay, '+');
    if (!truncated && opts.bothStrands) scan(revComp(hay), '-');
    out.sort(function (a, b) { return a.pos - b.pos || (a.strand < b.strand ? -1 : 1); });
    Object.defineProperty(out, 'truncated', { value: truncated, enumerable: false });
    return out;
  }

  /* =========================================================
   * 15. 双酶切兼容性速查
   * ========================================================= */
  /**
   * 两种酶对同一模板的联合分析：切点数、末端形态、是否同尾兼容、联合片段。
   * @param enzA/enzB 酶对象（MolEnzymes 中的条目）
   */
  function doubleDigestSummary(seq, enzA, enzB, circular) {
    if (!enzA || !enzB) return null;
    var hitsA = findSites(seq, enzA, circular);
    var hitsB = findSites(seq, enzB, circular);
    var endsA = enzymeEnds(enzA);
    var endsB = enzymeEnds(enzB);
    var d = digest(seq, [enzA, enzB], circular);
    var notes = [];
    if (enzA.name === enzB.name) notes.push('两种酶相同，等同单酶切。');
    var compatible = endsA && endsB && !endsA.degenerate && !endsB.degenerate &&
      endsA.type !== 'blunt' && endsB.type !== 'blunt' ? endsCompatible(endsA, endsB) :
      (endsA && endsB && endsA.type === 'blunt' && endsB.type === 'blunt');
    if (hitsA.length === 0) notes.push(enzA.name + ' 在模板内无切点。');
    if (hitsB.length === 0) notes.push(enzB.name + ' 在模板内无切点。');
    if (hitsA.length > 3) notes.push(enzA.name + ' 切点较多（' + hitsA.length + ' 个），片段会偏碎。');
    if (hitsB.length > 3) notes.push(enzB.name + ' 切点较多（' + hitsB.length + ' 个），片段会偏碎。');
    if (compatible && enzA.name !== enzB.name) notes.push('两者产生互补黏端（同尾/同裂）：切开后可互相连接，但连接点两种酶都再也切不动，适合「消灭位点」策略，不适合定向克隆。');
    if (endsA && endsB && endsA.type !== 'blunt' && endsB.type !== 'blunt' && !compatible && enzA.name !== enzB.name)
      notes.push('末端互不互补：可定向克隆，插入方向唯一。');
    if ((endsA && endsA.type === 'blunt') || (endsB && endsB.type === 'blunt'))
      notes.push('含平末端：连接效率低，方向不唯一，建议提高插入比例与连接时间。');
    return {
      cutsA: hitsA.length, cutsB: hitsB.length,
      endsA: endsA, endsB: endsB, compatible: !!compatible,
      fragments: d.fragments, notes: notes
    };
  }

  function formatSeq(seq, block) {
    block = block || 60;
    var out = [];
    for (var i = 0; i < seq.length; i += block) {
      out.push(seq.substr(i, block));
    }
    return out.join('\n');
  }

  function randomSeq(len, gc) {
    gc = gc === undefined ? 50 : gc;
    var s = '';
    for (var i = 0; i < len; i++) {
      var r = Math.random();
      if (r < gc / 100) s += (Math.random() < 0.5 ? 'G' : 'C');
      else s += (Math.random() < 0.5 ? 'A' : 'T');
    }
    return s;
  }

  function fmt(n, digits) {
    if (n === null || n === undefined || isNaN(n)) return '—';
    digits = digits === undefined ? 2 : digits;
    var v = Number(n);
    if (Math.abs(v) >= 1e6) return v.toExponential(3);
    return v.toFixed(digits).replace(/\.0+$/, '').replace(/(\.\d*?)0+$/, '$1');
  }

  global.MolCore = {
    WATER: WATER,
    CODON_TABLE: CODON_TABLE,
    AA_RESIDUE: AA_RESIDUE,
    DNA_RESIDUE: DNA_RESIDUE,
    RNA_RESIDUE: RNA_RESIDUE,
    IUPAC: IUPAC,
    COMPLEMENT: COMPLEMENT,
    OD_FACTOR: OD_FACTOR,
    parseFasta: parseFasta,
    parseGenBankFeatures: parseGenBankFeatures,
    extractFeature: extractFeature,
    parseLocation: parseLocation,
    detectFormat: detectFormat,
    detectType: detectType,
    isValidNucleic: isValidNucleic,
    baseCount: baseCount,
    gcContent: gcContent,
    gcWindows: gcWindows,
    revComp: revComp,
    transcribe: transcribe,
    reverseTranscribe: reverseTranscribe,
    translate: translate,
    findORFs: findORFs,
    findORFsCircular: findORFsCircular,
    tmWallace: tmWallace,
    tmGC: tmGC,
    tmNearestNeighbor: tmNearestNeighbor,
    nnSum: nnSum,
    nnDeltaG37: nnDeltaG37,
    hairpinDeltaG: hairpinDeltaG,
    dimerDeltaG: dimerDeltaG,
    dimerTouches3Prime: dimerTouches3Prime,
    mwNucleic: mwNucleic,
    mwProtein: mwProtein,
    fromOD: fromOD,
    concNgUlToNm: concNgUlToNm,
    concNmToNgUl: concNmToNgUl,
    findSites: findSites,
    fragmentsFromCuts: fragmentsFromCuts,
    digest: digest,
    gcClamp: gcClamp,
    pcrProducts: pcrProducts,
    pcrProductSeq: pcrProductSeq,
    enzymeEnds: enzymeEnds,
    endsCompatible: endsCompatible,
    digestOneDetailed: digestOneDetailed,
    ligateInsert: ligateInsert,
    designPrimers: designPrimers,
    alignPair: alignPair,
    alignMultiple: alignMultiple,
    proteinProperties: proteinProperties,
    proteinNetCharge: proteinNetCharge,
    reverseTranslate: reverseTranslate,
    ligationAmounts: ligationAmounts,
    designMutPrimers: designMutPrimers,
    findMotifs: findMotifs,
    doubleDigestSummary: doubleDigestSummary,
    formatSeq: formatSeq,
    randomSeq: randomSeq,
    fmt: fmt
  };
})(typeof window !== 'undefined' ? window : this);
