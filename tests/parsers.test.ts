/**
 * parsers.test.ts — parsers.ts 的 Node 24 原生测试（node:parsers.test.ts 直接运行）。
 *
 * 覆盖 qPCR-app/tests/fixtures/ 下全部文本类样例（xlsx 样例跳过：
 * Excel 解码在应用层用 SheetJS 完成，之后复用同一套网格函数 parseGrid）。
 *
 * 每个样例同时做两层验证：
 *   1. 硬编码行为断言（仪器识别、行数、NTC/Omit 过滤、未检出→null 等）；
 *   2. 与同目录 python_baseline.json（由 make_baseline.py 生成）逐字段对比，
 *      保证与 Python 版 parsers.py 零行为差异。
 *
 * 运行：node parsers.test.ts
 */
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

import {
  InstrumentParseError,
  SYNONYMS,
  CT_FALLBACK,
  UNDETERMINED,
  BINARY_HINTS,
  decodeText,
  textToGrid,
  chooseSep,
  parseDelimited,
  looksLikeHtml,
  htmlToGrid,
  scoreRow,
  detectHeaderRow,
  classifyInstrument,
  mapColumns,
  cleanCt,
  dedupHeader,
  gridToRawTable,
  lc480Channel,
  parseGrid,
  parseTextFile,
  parseTextFiles,
  parseManualGrid,
  autoGroupSamples,
  targetsInAllSamples,
  pivotToWide,
} from "../src/lib/parsers.ts";

const HERE = dirname(fileURLToPath(import.meta.url));
const FIXTURES = join(HERE, "..", "..", "qPCR-app", "tests", "fixtures");

let passed = 0;
function ok(name) {
  passed += 1;
  console.log(`✅ ${passed}. ${name}`);
}

function readFixture(name) {
  return readFileSync(join(FIXTURES, name));
}

function readFixtureText(name) {
  return decodeText(new Uint8Array(readFixture(name)));
}

const baseline = JSON.parse(readFileSync(join(HERE, "python_baseline.json"), "utf-8"));

function approxEqual(a, b, eps = 1e-9) {
  if (a === null && b === null) return true;
  if (a === null || b === null) return false;
  return Math.abs(a - b) < eps;
}

/** rows 与基准逐字段对比（ct 浮点容差 1e-9）。 */
function assertRowsMatch(actual, expected, label) {
  assert.equal(actual.length, expected.length, `${label}: 行数不一致`);
  for (let i = 0; i < actual.length; i++) {
    const a = actual[i];
    const e = expected[i];
    assert.equal(a.well, e.well, `${label}[${i}].well`);
    assert.equal(a.sample, e.sample, `${label}[${i}].sample`);
    assert.equal(a.target, e.target, `${label}[${i}].target`);
    assert.ok(approxEqual(a.ct, e.ct), `${label}[${i}].ct: ${a.ct} != ${e.ct}`);
    assert.equal(a.source_file, e.source_file, `${label}[${i}].source_file`);
  }
}

function assertReportMatch(actual, expected, label) {
  assert.deepEqual(actual, expected, `${label}: report 与 Python 基准不一致`);
}

// ==========================================================================
// 0. 常量与 Python 版逐字符一致（含西里尔 "cт"、全角破折号 "—" 等陷阱）
// ==========================================================================
assert.deepEqual(SYNONYMS, baseline.constants.SYNONYMS, "SYNONYMS 与 Python 版不一致");
assert.deepEqual(CT_FALLBACK, baseline.constants.CT_FALLBACK);
assert.deepEqual([...UNDETERMINED].sort(), baseline.constants.UNDETERMINED);
assert.deepEqual(BINARY_HINTS, baseline.constants.BINARY_HINTS);
ok("常量表（同义词/兜底/未检出/二进制提示）与 Python 版逐字符一致");

// ==========================================================================
// 1. ABI v1.x：Tab + 元信息块 + Results 节 + Undetermined/NTC/Omit
// ==========================================================================
{
  const text = readFixtureText("abi_v1_quantstudio.txt");
  const { rows, report } = parseTextFile("abi_v1_quantstudio.txt", text);
  assert.equal(report.status, "ok");
  assert.equal(report.instrument, "ABI/Thermo QuantStudio");
  assert.equal(rows.length, 12, "NTC 与 Omit=TRUE 行应被过滤");
  const a4 = rows.find((r) => r.well === "4");
  assert.equal(a4.ct, null, "Undetermined 应为 null");
  assert.deepEqual(report.targets, ["GAPDH", "MYC"]);
  const w1 = rows.find((r) => r.well === "1");
  assert.ok(approxEqual(w1.ct, 12.118));
  assert.equal(report.needs_target_name, false);
  assertRowsMatch(rows, baseline.files["abi_v1_quantstudio.txt"].rows, "abi_v1");
  assertReportMatch(report, baseline.files["abi_v1_quantstudio.txt"].reports[0], "abi_v1");
  // 表头行不写死：detectHeaderRow 应定位到第 16 行（0 基 15）
  const grid = textToGrid(text);
  assert.equal(detectHeaderRow(grid).headerIdx, 15);
}
ok("ABI v1.x（表头扫描定位 / NTC / Omit / Undetermined）解析通过");

// ==========================================================================
// 1b. ABI v2.x 风格导出（Well Position/Sample/Target/Cq）
// ==========================================================================
{
  const text = readFixtureText("abi_v2_quantstudio.txt");
  const { rows, report } = parseTextFile("abi_v2_quantstudio.txt", text);
  assert.equal(report.status, "ok");
  assert.equal(rows.length, 6);
  assert.ok(approxEqual(rows[1].ct, 23.2));
  assertRowsMatch(rows, baseline.files["abi_v2_quantstudio.txt"].rows, "abi_v2");
  assertReportMatch(report, baseline.files["abi_v2_quantstudio.txt"].reports[0], "abi_v2");
}
ok("ABI v2.x（Sample/Target/Cq 列名）解析通过");

// ==========================================================================
// 2. Bio-Rad CFX：分号分隔 + 欧区小数逗号 + 首列多余分隔符 + N/A
// ==========================================================================
{
  const text = readFixtureText("biorad_cfx.csv");
  const grid = textToGrid(text);
  // 分隔符嗅探应选中分号
  assert.equal(chooseSep(text.replace(/\r\n/g, "\n").split("\n")), ";");
  const { rows, report } = parseTextFile("biorad_cfx.csv", text);
  assert.equal(report.instrument, "Bio-Rad CFX");
  assert.equal(rows.length, 5, "Content=NTC 行应被过滤");
  const a02 = rows.find((r) => r.well === "A02");
  assert.ok(approxEqual(a02.ct, 23.57), "小数逗号 23,57 应转为 23.57");
  const a03 = rows.find((r) => r.well === "A03");
  assert.equal(a03.ct, null, "N/A 应为 null");
  // 首列多余分隔符产生的空 col_0 应被丢弃
  const { headerIdx } = detectHeaderRow(grid);
  const raw = gridToRawTable(grid, headerIdx);
  assert.equal(raw.columns[0], "Well", "空首列应被移除");
  assertRowsMatch(rows, baseline.files["biorad_cfx.csv"].rows, "biorad");
  assertReportMatch(report, baseline.files["biorad_cfx.csv"].reports[0], "biorad");
}
ok("Bio-Rad CFX（分号 + 小数逗号 + 空首列 + NTC + N/A）解析通过");

// ==========================================================================
// 3. Roche LightCycler 480：元信息行 Filter Combination → 通道名
// ==========================================================================
{
  for (const [fname, channel] of [
    ["lc480_GAPDH_533-610.txt", "533-610"],
    ["lc480_MYC_465-510.txt", "465-510"],
  ]) {
    const text = readFixtureText(fname);
    const grid = textToGrid(text);
    const { headerIdx } = detectHeaderRow(grid);
    assert.equal(lc480Channel(grid, headerIdx), channel);
    const { rows, report } = parseTextFile(fname, text);
    assert.equal(report.instrument, "Roche LightCycler 480");
    assert.equal(rows.length, 5);
    assert.equal(report.needs_target_name, true, "LC480 单通道文件应提示命名靶标");
    assert.equal(report.default_target, `LC480通道${channel}`);
    assert.deepEqual(report.targets, [`LC480通道${channel}`]);
    const ntc = rows.find((r) => r.sample === "NTC");
    assert.equal(ntc.ct, null, "空 Cp 应为 null");
    assertRowsMatch(rows, baseline.files[fname].rows, fname);
    assertReportMatch(report, baseline.files[fname].reports[0], fname);
  }
}
ok("LightCycler 480（Filter Combination 通道 + needsTargetName + 空 Cp）解析通过");

// ==========================================================================
// 4. LightCycler 96：Gene 列区分多基因
// ==========================================================================
{
  const text = readFixtureText("lc96.txt");
  const { rows, report } = parseTextFile("lc96.txt", text);
  assert.equal(report.instrument, "Roche LightCycler 96");
  assert.equal(rows.length, 6);
  assert.deepEqual(report.targets, ["GAPDH", "MYC"]);
  assertRowsMatch(rows, baseline.files["lc96.txt"].rows, "lc96");
  assertReportMatch(report, baseline.files["lc96.txt"].reports[0], "lc96");
}
ok("LightCycler 96（Gene 列）解析通过");

// ==========================================================================
// 5. Qiagen Rotor-Gene Q：数字孔位 + Type=NTC 过滤 + target 取文件名
// ==========================================================================
{
  const text = readFixtureText("rotor_gene_GAPDH.csv");
  const { rows, report } = parseTextFile("rotor_gene_GAPDH.csv", text);
  assert.equal(report.instrument, "Qiagen Rotor-Gene Q");
  assert.equal(rows.length, 4, "Type=NTC 行应被过滤");
  assert.equal(rows[0].target, "rotor_gene_GAPDH", "target 应取文件名（去扩展名）");
  assert.equal(report.needs_target_name, true);
  assert.ok(approxEqual(rows[0].ct, 14.21));
  assertRowsMatch(rows, baseline.files["rotor_gene_GAPDH.csv"].rows, "rotor");
  assertReportMatch(report, baseline.files["rotor_gene_GAPDH.csv"].reports[0], "rotor");
}
ok("Rotor-Gene Q（数字孔位 + NTC 过滤 + target 取文件名）解析通过");

// ==========================================================================
// 6. 国产仪器：中文表头 + GBK 编码（decodeText 回退 gb18030）+ 中文未检出
// ==========================================================================
{
  const bytes = new Uint8Array(readFixture("gentier_国产.csv"));
  // 该文件是 GBK 编码，utf-8 严格解码应失败并回退 gb18030
  assert.throws(() => new TextDecoder("utf-8", { fatal: true }).decode(bytes));
  const text = decodeText(bytes);
  assert.ok(text.includes("孔位"), "GBK 解码后应包含中文表头");
  const { rows, report } = parseTextFile("gentier_国产.csv", text);
  assert.equal(report.status, "ok");
  assert.equal(report.instrument, "国产仪器/通用格式（中文表头）");
  assert.equal(rows.length, 5);
  const s2 = rows.filter((r) => r.sample === "样本2");
  const undet = s2.find((r) => r.target === "目的基因");
  assert.equal(undet.ct, null, "中文「未检出」应为 null");
  assertRowsMatch(rows, baseline.files["gentier_国产.csv"].rows, "gentier");
  assertReportMatch(report, baseline.files["gentier_国产.csv"].reports[0], "gentier");
}
ok("国产仪器（中文表头 + GBK 回退 + 未检出）解析通过");

// ==========================================================================
// 7. HTML 伪 xls：looksLikeHtml 识别 + 表格提取（oldsoftware_results.xls）
// ==========================================================================
{
  const text = readFixtureText("oldsoftware_results.xls");
  assert.ok(looksLikeHtml(text), "应识别为 HTML 伪 xls");
  const grid = htmlToGrid("oldsoftware_results.xls", text);
  assert.deepEqual(grid[0], ["Well", "Sample Name", "Target", "CT"]);
  const { rows, report } = parseTextFile("oldsoftware_results.xls", text);
  assert.equal(report.status, "ok");
  assert.equal(rows.length, 4);
  const a4 = rows.find((r) => r.well === "A4");
  assert.equal(a4.ct, null, "Undetermined 应为 null");
  assertRowsMatch(rows, baseline.files["oldsoftware_results.xls"].rows, "oldsoftware");
  assertReportMatch(report, baseline.files["oldsoftware_results.xls"].reports[0], "oldsoftware");
}
ok("HTML 伪 xls（表格提取 + Undetermined）解析通过");

// ==========================================================================
// 8. 错误路径：二进制工程文件提示 / 无表头 / 表头后无数据 / 手动行号越界
// ==========================================================================
{
  assert.throws(
    () => parseTextFile("exp.eds", "whatever"),
    (e) => e instanceof InstrumentParseError && e.message.includes("Export")
  );
  assert.throws(
    () => parseTextFile("run.pcrd", "whatever"),
    (e) => e instanceof InstrumentParseError && e.message.includes("导出")
  );
  assert.throws(
    () => parseTextFile("random.txt", "hello\tworld\n1\t2\n"),
    (e) => e instanceof InstrumentParseError && e.message.includes("未能定位结果表头行")
  );
  assert.throws(
    () => parseTextFile("nodata.csv", "Well,Sample,Cq\n"),
    (e) => e instanceof InstrumentParseError && e.message.includes("表头之后没有数据行")
  );
  const grid = textToGrid(readFixtureText("abi_v1_quantstudio.txt"));
  assert.throws(
    () => parseManualGrid(grid, "x.txt", 999, { ct: "CT" }),
    (e) => e instanceof InstrumentParseError && e.message.includes("超出范围")
  );
  // 只有 Cq Mean 列时走 CT_FALLBACK 兜底
  const fb = mapColumns(["Well", "Sample", "Cq Mean"]);
  assert.equal(fb.ct, "Cq Mean", "无精确 Ct 列时应兜底均值列");
  // 缺 Ct 列与缺样本/孔位列的报错
  assert.throws(
    () => parseGrid([["Well", "Sample"], ["A1", "S1"]], "x.csv"),
    (e) => e instanceof InstrumentParseError && e.message.includes("未能定位结果表头行")
  );
  // parseTextFiles 单文件失败不中断其他文件
  const multi = parseTextFiles([
    { name: "bad.eds", text: "binary" },
    { name: "lc96.txt", text: readFixtureText("lc96.txt") },
  ]);
  assert.equal(multi.reports[0].status, "error");
  assert.equal(multi.reports[1].status, "ok");
  assert.equal(multi.rows.length, 6);
}
ok("错误路径（.eds/.pcrd 提示、无表头、无数据、越界、兜底、失败不中断）通过");

// ==========================================================================
// 9. 多文件合并：全部文本样例（xlsx 样例跳过）与 Python 基准一致
// ==========================================================================
{
  const names = [
    "abi_v1_quantstudio.txt",
    "abi_v2_quantstudio.txt",
    "biorad_cfx.csv",
    "lc480_GAPDH_533-610.txt",
    "lc480_MYC_465-510.txt",
    "lc96.txt",
    "rotor_gene_GAPDH.csv",
    "gentier_国产.csv",
    "oldsoftware_results.xls",
  ];
  const files = names.map((n) => ({ name: n, text: readFixtureText(n) }));
  const { rows, reports } = parseTextFiles(files);
  assert.equal(reports.length, names.length);
  assert.ok(reports.every((r) => r.status === "ok"));
  assertRowsMatch(rows, baseline.multi.rows, "multi");
  assert.deepEqual(reports, baseline.multi.reports, "multi reports 与基准不一致");
  // 合并结构：source_file 保留来源
  const sources = new Set(rows.map((r) => r.source_file));
  assert.equal(sources.size, names.length);
}
ok("多文件合并（9 个文本样例，52 行长表 + 逐文件报告）与 Python 基准一致");

// ==========================================================================
// 10. 手动模式：指定表头行 + 列名
// ==========================================================================
{
  const text = readFixtureText("abi_v1_quantstudio.txt");
  const grid = textToGrid(text);
  const { rows, report } = parseManualGrid(grid, "abi_v1_quantstudio.txt", 15, {
    well: "Well",
    sample: "Sample Name",
    target: "Target Name",
    ct: "CT",
  });
  assert.equal(report.status, "ok");
  assert.equal(rows.length, 12);
  assert.ok(report.instrument.includes("手动指定"));
  assertRowsMatch(rows, baseline.manual.rows, "manual");
  assert.deepEqual(report, baseline.manual.reports[0], "manual report 与基准不一致");
}
ok("手动模式（指定表头行 + 列映射）与 Python 基准一致");

// ==========================================================================
// 11. 下游辅助：自动分组 / 管家基因候选 / 长表转宽表
// ==========================================================================
{
  const names = [
    "abi_v1_quantstudio.txt", "abi_v2_quantstudio.txt", "biorad_cfx.csv",
    "lc480_GAPDH_533-610.txt", "lc480_MYC_465-510.txt", "lc96.txt",
    "rotor_gene_GAPDH.csv", "gentier_国产.csv", "oldsoftware_results.xls",
  ];
  const { rows } = parseTextFiles(names.map((n) => ({ name: n, text: readFixtureText(n) })));
  const samples = [...new Set(rows.map((r) => r.sample))].sort();

  assert.deepEqual(autoGroupSamples(samples, "identity"), baseline.helpers.group_identity);
  assert.deepEqual(autoGroupSamples(samples, "strip_suffix"), baseline.helpers.group_strip);
  assert.deepEqual(autoGroupSamples(samples, "prefix_sep"), baseline.helpers.group_prefix);
  assert.deepEqual(
    autoGroupSamples(["NC_1", "NC_2", "TREAT_1", "TREAT_2"], "strip_suffix"),
    { NC_1: "NC", NC_2: "NC", TREAT_1: "TREAT", TREAT_2: "TREAT" }
  );

  const hk = targetsInAllSamples(rows);
  assert.deepEqual(hk, baseline.helpers.hk_targets, "管家基因候选排序与基准不一致");
  assert.equal(hk[0], "GAPDH", "覆盖全部样本的 GAPDH 应置顶");

  const wide = pivotToWide(rows);
  assert.deepEqual(wide.columns, baseline.helpers.wide_columns);
  assert.equal(wide.rows.length, baseline.helpers.wide_rows.length);
  for (let i = 0; i < wide.rows.length; i++) {
    const a = wide.rows[i];
    const e = baseline.helpers.wide_rows[i];
    for (const col of wide.columns) {
      if (col === "sample") {
        assert.equal(a[col], e[col], `wide[${i}].sample`);
      } else {
        assert.ok(
          approxEqual(a[col] ?? null, e[col] ?? null),
          `wide[${i}].${col}: ${a[col]} != ${e[col]}`
        );
      }
    }
  }
}
ok("自动分组 / 管家基因候选 / 长表转宽表 与 Python 基准一致");

// ==========================================================================
// 12. 底层小工具单元行为
// ==========================================================================
{
  // cleanCt：欧区小数逗号、未检出集合、非法值
  assert.ok(approxEqual(cleanCt("23,5729"), 23.5729));
  assert.ok(approxEqual(cleanCt(" 12.1 "), 12.1));
  assert.equal(cleanCt("Undetermined"), null);
  assert.equal(cleanCt("未检出"), null);
  assert.equal(cleanCt("—"), null);
  assert.equal(cleanCt("abc"), null);
  assert.equal(cleanCt("1234,5"), null, "整数部分超过 3 位不当作小数逗号");
  // dedupHeader：空名补 col_j，重名加后缀
  assert.deepEqual(dedupHeader(["Well", "", "Well", "CT"]), ["Well", "col_1", "Well_1", "CT"]);
  // parseDelimited：引号字段与转义引号
  assert.deepEqual(parseDelimited('a,"b,c","d""e"\n1,2,3', ","), [
    ["a", "b,c", 'd"e'],
    ["1", "2", "3"],
  ]);
  // scoreRow / classifyInstrument 规则链
  const { score, fields } = scoreRow(["Well", "Sample Name", "Target Name", "CT"]);
  assert.equal(score, 4);
  assert.ok(fields.has("ct") && fields.has("well") && fields.has("sample") && fields.has("target"));
  assert.equal(classifyInstrument(["Well", "Target", "Ct"]), "ABI/Thermo QuantStudio (v2.x)");
  assert.equal(classifyInstrument(["孔位", "样本名称", "Ct值"]), "国产仪器/通用格式（中文表头）");
  assert.equal(classifyInstrument(["Well", "Foo", "Ct"]), "通用格式（按列名同义词识别）");
  // decodeText：utf-8 BOM 正常处理
  const bom = new Uint8Array([0xef, 0xbb, 0xbf, 0x61, 0x62]);
  assert.equal(decodeText(bom), "ab");
}
ok("底层工具（cleanCt/dedupHeader/CSV 引号/规则链/BOM）单元行为通过");

console.log(`\n全部 ${passed} 组测试通过 🎉（所有文本 fixtures 与 Python 基准零差异）`);
