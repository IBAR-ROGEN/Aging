/**
 * Build the Vite app, serve preview, screenshot #dashboard-figure-root to analysis/.
 */
import { spawn } from "node:child_process";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { chromium } from "playwright";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const appRoot = path.resolve(__dirname, "..");
const repoRoot = path.resolve(appRoot, "..", "..");
const outDir = path.join(repoRoot, "analysis");
const outFile = path.join(outDir, "dashboard_figure_mockup.png");

function run(cmd, args, opts) {
  return new Promise((resolve, reject) => {
    const p = spawn(cmd, args, { ...opts, stdio: "inherit" });
    p.on("error", reject);
    p.on("close", (code) =>
      code === 0 ? resolve() : reject(new Error(`${cmd} exited ${code}`)),
    );
  });
}

async function waitForHttp(url, timeoutMs = 60_000) {
  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {
    try {
      const res = await fetch(url);
      if (res.ok) {
        return;
      }
    } catch {
      /* retry */
    }
    await new Promise((r) => setTimeout(r, 250));
  }
  throw new Error(`Timeout waiting for ${url}`);
}

fs.mkdirSync(outDir, { recursive: true });

await run("npm", ["run", "build"], { cwd: appRoot, shell: true });

const preview = spawn("npm", ["run", "preview:static"], {
  cwd: appRoot,
  stdio: ["ignore", "pipe", "pipe"],
  shell: true,
});

try {
  await waitForHttp("http://127.0.0.1:4173/");
  const browser = await chromium.launch();
  const page = await browser.newPage({
    viewport: { width: 1600, height: 980 },
    deviceScaleFactor: 2,
  });
  await page.goto("http://127.0.0.1:4173/", { waitUntil: "networkidle" });
  await page.locator("#dashboard-figure-root").waitFor({ state: "visible" });
  await new Promise((r) => setTimeout(r, 2000));
  await page.locator("#dashboard-figure-root").screenshot({
    path: outFile,
    type: "png",
  });
  await browser.close();
  console.log(`Wrote ${outFile}`);
} finally {
  preview.kill("SIGTERM");
  await new Promise((r) => setTimeout(r, 500));
}
