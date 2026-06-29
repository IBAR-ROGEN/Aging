# Dashboard figure render (React → PNG)

Vite + Tailwind + Recharts build of [`../DashboardFigureMockup.tsx`](../DashboardFigureMockup.tsx) for a **pixel-accurate** manuscript screenshot.

```bash
npm install
npx playwright install chromium
npm run capture
```

Writes `../../figures/dashboard_figure_mockup.png` (same path as the matplotlib version from `scripts/figures/render_dashboard_figure_mockup.py`).

```bash
npm run dev
```

Preview at `http://localhost:5173` during design iterations.
