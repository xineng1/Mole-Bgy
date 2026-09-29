#!/usr/bin/env node
/*
 * MolBench → GitHub 推送脚本（走 api.github.com，绕开被拦截的 github.com 主站）
 *
 * 用法：
 *   GITHUB_TOKEN=ghp_xxx node tools/push_to_github.mjs xineng1/Mole-Bgy
 *   node tools/push_to_github.mjs xineng1/Mole-Bgy --token ghp_xxx
 *   node tools/push_to_github.mjs xineng1/Mole-Bgy --dry-run     # 只打印将要上传的文件，不发请求
 *
 * Token 要求（二选一）：
 *   · Fine-grained token：Repository access 选中该仓库，权限 Contents = Read and write
 *   · Classic token：勾选 repo 范围
 *
 * 实现：Git Data API（blobs → tree → commit → ref），空仓库也能一次推成，
 *       不依赖本地 git 凭据，也不受 github.com:443 被拦影响（api.github.com 可直连）。
 */
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(__dirname, '..');
const API = 'https://api.github.com';

// ---------- 参数 ----------
const argv = process.argv.slice(2);
const flags = argv.filter(a => a.startsWith('--'));
const positional = argv.filter(a => !a.startsWith('--'));
const DRY = flags.includes('--dry-run');
const tokenFlagIdx = argv.indexOf('--token');
const TOKEN = process.env.GITHUB_TOKEN || process.env.GH_TOKEN ||
  (tokenFlagIdx >= 0 ? argv[tokenFlagIdx + 1] : '') || '';
const REPO = positional[0] || 'xineng1/Mole-Bgy';
const BRANCH = 'main';

// ---------- 忽略规则 ----------
const IGNORE_DIRS = new Set(['.git', 'node_modules', 'dist', 'build', '.tmp', 'tmp']);
const IGNORE_FILES = new Set(['Thumbs.db', '.DS_Store', 'desktop.ini']);
const IGNORE_EXT = /\.(log|tmp|swp)$/i;
const IGNORE_NAME = /(~|\.bak)$/;

function collect(dir, base = '') {
  const out = [];
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const rel = base ? base + '/' + entry.name : entry.name;
    if (entry.isDirectory()) {
      if (IGNORE_DIRS.has(entry.name)) continue;
      out.push(...collect(path.join(dir, entry.name), rel));
    } else {
      if (IGNORE_FILES.has(entry.name) || IGNORE_EXT.test(entry.name) || IGNORE_NAME.test(entry.name)) continue;
      out.push(rel);
    }
  }
  return out.sort();
}

// ---------- GitHub API ----------
async function api(method, url, body) {
  const res = await fetch(API + url, {
    method,
    headers: {
      'Authorization': 'Bearer ' + TOKEN,
      'Accept': 'application/vnd.github+json',
      'X-GitHub-Api-Version': '2022-11-28',
      'User-Agent': 'MolBench-push-script',
      ...(body ? { 'Content-Type': 'application/json' } : {})
    },
    body: body ? JSON.stringify(body) : undefined
  });
  const text = await res.text();
  let json = null;
  try { json = text ? JSON.parse(text) : null; } catch { /* 非 JSON 响应 */ }
  if (!res.ok) {
    const msg = (json && (json.message || json.error)) || text.slice(0, 200);
    const err = new Error(`${method} ${url} → HTTP ${res.status}: ${msg}`);
    err.status = res.status;
    err.body = json;
    throw err;
  }
  return json;
}

const COMMIT_MESSAGE = `MolBench v2：离线分子生物实验台

零依赖单页工具，双击 index.html 即可离线使用。9 个功能面板：

- 序列概览：GC 滑动窗口曲线、序列条形码、分子量与 OD260 定量
- 序列操作 / ORF：反向互补、六框翻译、ORF 查找（含嵌套 ORF、环状跨接缝、负链坐标映射）
- 蛋白 & 密码子：pI / GRAVY / 不稳定指数 / 脂肪族指数 / 消光系数；E. coli 偏好密码子反向翻译
- 酶切分析 & 虚拟凝胶：51 种酶、简并位点、环状跨接缝切点、质粒环形图谱、SC/LIN/OC 构象带、
  Motif（IUPAC）搜索、双酶切兼容性速查
- 引物 & PCR：自动设计（Tm/GC/ΔG 二级结构/模板唯一性打分）、NN Tm + 盐校正、Q5 型定点突变引物
- 序列比对：Needleman-Wunsch 双序列 + 中心星多序列比对
- 虚拟克隆：黏端形态推导与同尾酶判定、连接反应摩尔比用量计算
- 实验计算器：稀释、质量⇄摩尔、PCR 体系配制、rpm⇄×g、菌浓估算
- 说明 / 算法出处

自检：8 套脚本约 400 项断言（check / dom_smoke / feature_audit / bug_hunt ×2 / cross_validate / lint / profile）
+ 真实 Chrome 像素自检，全部通过。`;

async function main() {
  const files = collect(ROOT);
  const totalBytes = files.reduce((s, f) => s + fs.statSync(path.join(ROOT, f)).size, 0);

  console.log(`仓库    : ${REPO}（分支 ${BRANCH}）`);
  console.log(`文件数  : ${files.length} 个，合计 ${(totalBytes / 1024).toFixed(1)} KB`);
  console.log(`Token   : ${TOKEN ? '已提供（' + TOKEN.slice(0, 7) + '…）' : '未提供'}`);
  console.log('');
  if (DRY) {
    console.log('[dry-run] 将上传以下文件：');
    files.forEach(f => console.log('  ' + f));
    console.log('\n[dry-run] 未发送任何请求。去掉 --dry-run 即真正推送。');
    return;
  }
  if (!TOKEN) {
    console.error('缺少 token。请设置环境变量 GITHUB_TOKEN，或用 --token 传入。');
    process.exit(2);
  }

  // 1. 目标分支当前状态
  let parentSha = null;
  try {
    const ref = await api('GET', `/repos/${REPO}/git/ref/heads/${BRANCH}`);
    parentSha = ref.object.sha;
    console.log(`分支已存在，父提交 ${parentSha.slice(0, 8)}`);
  } catch (e) {
    if (e.status === 404) console.log('分支尚不存在（空仓库），将创建首个提交');
    else throw e;
  }

  // 2. 逐个创建 blob
  const tree = [];
  let done = 0;
  for (const f of files) {
    const content = fs.readFileSync(path.join(ROOT, f));
    const blob = await api('POST', `/repos/${REPO}/git/blobs`, {
      content: content.toString('base64'),
      encoding: 'base64'
    });
    tree.push({ path: f, mode: '100644', type: 'blob', sha: blob.sha });
    done++;
    if (done % 5 === 0 || done === files.length) console.log(`  已上传 blob ${done}/${files.length}`);
  }

  // 3. 建树
  const treeRes = await api('POST', `/repos/${REPO}/git/trees`, { tree });
  console.log('树对象', treeRes.sha.slice(0, 8));

  // 4. 提交
  const commit = await api('POST', `/repos/${REPO}/git/commits`, {
    message: COMMIT_MESSAGE,
    tree: treeRes.sha,
    parents: parentSha ? [parentSha] : []
  });
  console.log('提交', commit.sha.slice(0, 8));

  // 5. 更新（或创建）分支引用
  if (parentSha) {
    await api('PATCH', `/repos/${REPO}/git/refs/heads/${BRANCH}`, { sha: commit.sha, force: false });
  } else {
    await api('POST', `/repos/${REPO}/git/refs`, { ref: `refs/heads/${BRANCH}`, sha: commit.sha });
  }

  const html = `https://github.com/${REPO}`;
  console.log('\n推送完成 ✔');
  console.log('提交页 : ' + html + '/commit/' + commit.sha);
  console.log('仓库页 : ' + html);
  console.log(`（共 ${files.length} 个文件，${(totalBytes / 1024).toFixed(1)} KB）`);
}

main().catch(e => {
  console.error('\n推送失败：' + e.message);
  if (e.status === 401) console.error('→ token 无效或已过期');
  if (e.status === 403) console.error('→ token 权限不足（需要 Contents: Read and write）或触发了速率限制');
  if (e.status === 404) console.error('→ 仓库不存在，或 token 未被授权访问该仓库');
  process.exit(1);
});
