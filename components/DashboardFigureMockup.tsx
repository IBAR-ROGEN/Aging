/**
 * Static high-fidelity mockup of the ROGEN Streamlit EDA dashboard (Epigenetic Clock tab).
 * Intended for one-off render + screenshot for manuscript figures.
 *
 * Dependencies: react, recharts, lucide-react, tailwindcss
 *
 * Example wrapper (ensure Inter is loaded, e.g. next/font or a link tag):
 *   <div className="min-h-screen bg-slate-100 p-8 font-sans antialiased">
 *     <DashboardFigureMockup />
 *   </div>
 */

import type { ReactNode } from "react";
import {
  Activity,
  Dna,
  FileSpreadsheet,
  FlaskConical,
  Settings2,
  SlidersHorizontal,
} from "lucide-react";
import {
  CartesianGrid,
  ComposedChart,
  Line,
  ResponsiveContainer,
  Scatter,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";

/** Deterministic ~50 points: y ≈ x with biological-looking spread (no useState). */
const MOCK_SCATTER_DATA: { x: number; y: number }[] = (() => {
  const out: { x: number; y: number }[] = [];
  let s = 42_069;
  for (let i = 0; i < 50; i++) {
    s = (s * 1_103_515_245 + 12_345) >>> 0;
    const u = s / 4_294_967_295;
    const t = i / 49;
    const x = 40 + 60 * t + (u - 0.5) * 3.5;
    const noise = Math.sin(i * 1.7 + u * 6) * 3.2 + (u - 0.48) * 4;
    const y = x + noise;
    out.push({
      x: Math.round(x * 100) / 100,
      y: Math.round(Math.min(102, Math.max(38, y)) * 100) / 100,
    });
  }
  return out;
})();

const IDENTITY_LINE: { x: number; y: number }[] = [
  { x: 40, y: 40 },
  { x: 100, y: 100 },
];

function MetricCard({
  label,
  value,
  sub,
}: {
  label: string;
  value: string;
  sub?: string;
}) {
  return (
    <div className="flex-1 rounded-lg border border-slate-200 bg-white px-5 py-4 shadow-sm">
      <p className="text-xs font-medium uppercase tracking-wide text-slate-500">
        {label}
      </p>
      <p className="mt-1 font-semibold tabular-nums text-slate-900 text-2xl tracking-tight">
        {value}
      </p>
      {sub ? (
        <p className="mt-0.5 text-xs text-slate-400">{sub}</p>
      ) : null}
    </div>
  );
}

function SidebarSectionTitle({ children }: { children: ReactNode }) {
  return (
    <h3 className="mb-2 flex items-center gap-2 text-xs font-semibold uppercase tracking-wider text-slate-500">
      <SlidersHorizontal className="h-3.5 w-3.5 text-slate-400" aria-hidden />
      {children}
    </h3>
  );
}

/** Visual only: dual-thumb range track for manuscript figure. */
function DualRangeSliderMock() {
  return (
    <div className="space-y-2">
      <div className="flex justify-between text-xs text-slate-600">
        <span className="font-medium text-slate-700">Age range</span>
        <span className="tabular-nums text-slate-500">40 — 100 years</span>
      </div>
      <div className="relative h-2 rounded-full bg-slate-200">
        <div
          className="absolute inset-y-0 rounded-full bg-sky-200/90"
          style={{ left: "4%", right: "8%" }}
        />
        <div
          className="absolute top-1/2 h-4 w-4 -translate-x-1/2 -translate-y-1/2 rounded-full border-2 border-sky-600 bg-white shadow-sm"
          style={{ left: "4%" }}
        />
        <div
          className="absolute top-1/2 h-4 w-4 -translate-x-1/2 -translate-y-1/2 rounded-full border-2 border-sky-600 bg-white shadow-sm"
          style={{ left: "92%" }}
        />
      </div>
    </div>
  );
}

function MultiselectMock({
  label,
  chips,
}: {
  label: string;
  chips: string[];
}) {
  return (
    <div className="space-y-1.5">
      <span className="text-xs font-medium text-slate-600">{label}</span>
      <div className="flex min-h-[2.5rem] flex-wrap items-center gap-1.5 rounded-md border border-slate-200 bg-white px-2 py-2 shadow-inner">
        {chips.map((c) => (
          <span
            key={c}
            className="inline-flex items-center rounded border border-sky-200 bg-sky-50 px-2 py-0.5 text-xs font-medium text-sky-800"
          >
            {c}
          </span>
        ))}
        <span className="ml-auto text-slate-300">
          <svg
            className="h-4 w-4"
            fill="none"
            stroke="currentColor"
            viewBox="0 0 24 24"
            aria-hidden
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M19 9l-7 7-7-7"
            />
          </svg>
        </span>
      </div>
    </div>
  );
}

function ToggleSwitchMock() {
  return (
    <div className="flex items-center justify-between gap-3 rounded-lg border border-slate-200 bg-slate-50/80 px-3 py-2.5">
      <span className="text-sm font-medium text-slate-700">
        Use Synthetic Mock Cohort
      </span>
      <div
        className="relative h-7 w-12 shrink-0 rounded-full bg-emerald-500 shadow-inner"
        aria-hidden
      >
        <div className="absolute right-0.5 top-0.5 h-6 w-6 rounded-full bg-white shadow-md ring-1 ring-black/5" />
      </div>
    </div>
  );
}

export default function DashboardFigureMockup() {
  return (
    <div
      id="dashboard-figure-root"
      className="mx-auto max-w-7xl overflow-hidden rounded-xl border border-slate-300/80 bg-white font-sans text-slate-800 shadow-xl antialiased"
      style={{ fontFamily: "Inter, ui-sans-serif, system-ui, sans-serif" }}
    >
      <div className="flex min-h-[640px]">
        {/* Sidebar ~20% */}
        <aside className="flex w-[20%] min-w-[220px] flex-col border-r border-slate-200 bg-slate-50/90">
          <div className="border-b border-slate-200 bg-white px-4 py-4">
            <div className="flex items-center gap-2 text-slate-800">
              <Settings2 className="h-5 w-5 text-sky-600" strokeWidth={2} />
              <h2 className="text-sm font-bold tracking-tight">
                ROGEN EDA Settings
              </h2>
            </div>
          </div>

          <div className="flex flex-1 flex-col gap-5 overflow-hidden px-4 py-4">
            <div className="space-y-2">
              <label className="flex items-center gap-1.5 text-xs font-medium text-slate-600">
                <FileSpreadsheet className="h-3.5 w-3.5 text-slate-400" />
                Merged cohort path
              </label>
              <div className="rounded-md border border-slate-200 bg-white px-3 py-2 font-mono text-xs text-slate-700 shadow-sm">
                data/merged_cohort.parquet
              </div>
            </div>

            <ToggleSwitchMock />

            <div className="h-px bg-slate-200" />

            <div>
              <SidebarSectionTitle>Global filters</SidebarSectionTitle>
              <div className="space-y-4">
                <DualRangeSliderMock />
                <MultiselectMock
                  label="Sex"
                  chips={["Female", "Male"]}
                />
                <MultiselectMock
                  label="Disease status"
                  chips={["Control", "Case", "Prodromal"]}
                />
              </div>
            </div>
          </div>
        </aside>

        {/* Main ~80% */}
        <main className="flex min-w-0 flex-1 flex-col bg-white">
          <header className="border-b border-slate-100 px-8 py-6">
            <h1 className="text-2xl font-bold tracking-tight text-slate-900 md:text-[1.65rem]">
              ROGEN Multi-Omics Exploratory Data Analysis
            </h1>
            <p className="mt-2 max-w-4xl text-sm leading-relaxed text-slate-600">
              Integrating whole-genome sequence-derived variants, DNA methylation–based
              age estimates, and structured clinical phenotypes to characterize
              aging-related biology in the merged cohort.
            </p>
          </header>

          {/* Tabs */}
          <div className="border-b border-slate-200 bg-slate-50/50 px-6">
            <nav
              className="flex gap-1"
              aria-label="Dashboard sections"
            >
              <button
                type="button"
                className="flex items-center gap-2 border-b-2 border-transparent px-4 py-3 text-sm font-medium text-slate-500"
                tabIndex={-1}
              >
                <FlaskConical className="h-4 w-4" strokeWidth={2} />
                Clinical &amp; Phenotypic
              </button>
              <button
                type="button"
                className="-mb-px flex items-center gap-2 border-b-2 border-sky-600 bg-white px-4 py-3 text-sm font-semibold text-sky-800 shadow-[0_-1px_0_0_white]"
                tabIndex={-1}
              >
                <Activity className="h-4 w-4 text-sky-600" strokeWidth={2} />
                Epigenetic Clock Validation
              </button>
              <button
                type="button"
                className="flex items-center gap-2 border-b-2 border-transparent px-4 py-3 text-sm font-medium text-slate-500"
                tabIndex={-1}
              >
                <Dna className="h-4 w-4" strokeWidth={2} />
                LA-SNPs
              </button>
            </nav>
          </div>

          <div className="flex flex-1 flex-col gap-6 px-8 py-6">
            <div>
              <h2 className="text-lg font-semibold text-slate-900">
                Epigenetic clock validation
              </h2>
              <p className="mt-1 text-sm text-slate-600">
                Chronological age versus predicted DNAm age with ordinary least squares
                fit; metrics summarize agreement for the filtered cohort.
              </p>
            </div>

            <div className="flex flex-wrap gap-4">
              <MetricCard
                label="Mean Absolute Error (MAE)"
                value="3.42 years"
              />
              <MetricCard
                label="Pearson Correlation (r)"
                value="0.89"
              />
              <MetricCard label="Cohort Size (N)" value="842" />
            </div>

            <div className="rounded-lg border border-slate-200 bg-white p-4 shadow-sm">
              <p className="mb-2 text-center text-sm font-medium text-slate-700">
                Chronological age vs. predicted DNAm age
              </p>
              <div className="h-96 w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <ComposedChart
                    margin={{ top: 8, right: 24, left: 8, bottom: 28 }}
                  >
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                      type="number"
                      dataKey="x"
                      domain={[40, 100]}
                      ticks={[40, 50, 60, 70, 80, 90, 100]}
                      tick={{ fill: "#475569", fontSize: 11 }}
                      label={{
                        value: "Chronological Age (years)",
                        position: "bottom",
                        offset: 12,
                        fill: "#334155",
                        fontSize: 12,
                        fontWeight: 500,
                      }}
                    />
                    <YAxis
                      type="number"
                      dataKey="y"
                      domain={[40, 100]}
                      ticks={[40, 50, 60, 70, 80, 90, 100]}
                      tick={{ fill: "#475569", fontSize: 11 }}
                      label={{
                        value: "Predicted DNAm Age (years)",
                        angle: -90,
                        position: "insideLeft",
                        fill: "#334155",
                        fontSize: 12,
                        fontWeight: 500,
                        dx: -6,
                      }}
                    />
                    <Tooltip
                      contentStyle={{
                        borderRadius: "6px",
                        border: "1px solid #e2e8f0",
                        fontSize: "12px",
                      }}
                      formatter={(v: number | string) => [
                        typeof v === "number" ? v.toFixed(2) : v,
                        "",
                      ]}
                      labelFormatter={() => "Sample"}
                    />
                    <Line
                      data={IDENTITY_LINE}
                      type="linear"
                      dataKey="y"
                      stroke="#64748b"
                      strokeWidth={2}
                      dot={false}
                      strokeDasharray="6 5"
                      isAnimationActive={false}
                      name="OLS (y = x)"
                    />
                    <Scatter
                      name="Samples"
                      data={MOCK_SCATTER_DATA}
                      fill="#0284c7"
                      fillOpacity={0.88}
                      isAnimationActive={false}
                    />
                  </ComposedChart>
                </ResponsiveContainer>
              </div>
              <p className="mt-2 text-center text-xs text-slate-500">
                Dashed line: identity / OLS reference (y = x). Points: independent
                samples (mock data for illustration).
              </p>
            </div>
          </div>
        </main>
      </div>
    </div>
  );
}
