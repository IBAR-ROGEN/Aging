# Dashboard figure render (React → PNG)

Vite + Tailwind + Recharts build of [`../DashboardFigureMockup.tsx`](../DashboardFigureMockup.tsx) for a **pixel-accurate** manuscript screenshot.

```bash
npm install
npx playwright install chromium
npm run capture
```

Writes `../../analysis/dashboard_figure_mockup.png` (overwrites the matplotlib version from `scripts/render_dashboard_figure_mockup.py` if you use the same path).

```bash
npm run dev
```

Preview at `http://localhost:5173` during design iterations.
