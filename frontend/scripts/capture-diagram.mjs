import { spawn } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium } from 'playwright';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const frontendRoot = path.resolve(__dirname, '..');
const repoRoot = path.resolve(frontendRoot, '..');
const outFile = path.join(repoRoot, 'figures', 'longevity_network_diagram.png');

function run(cmd, args, options) {
  return new Promise((resolve, reject) => {
    const child = spawn(cmd, args, {
      stdio: 'inherit',
      shell: process.platform === 'win32',
      ...options,
    });
    child.on('error', reject);
    child.on('close', (code) => {
      if (code === 0) {
        resolve();
      } else {
        reject(new Error(`${cmd} ${args.join(' ')} exited with ${code}`));
      }
    });
  });
}

function startPreview() {
  const child = spawn('npm', ['run', 'preview', '--', '--host', '127.0.0.1', '--port', '4173'], {
    cwd: frontendRoot,
    stdio: ['ignore', 'pipe', 'pipe'],
    shell: process.platform === 'win32',
  });
  return child;
}

async function waitForHttpOk(url, timeoutMs) {
  const start = Date.now();
  while (Date.now() - start < timeoutMs) {
    try {
      const res = await fetch(url, { signal: AbortSignal.timeout(2000) });
      if (res.ok) {
        return;
      }
    } catch {
      /* retry */
    }
    await new Promise((r) => setTimeout(r, 300));
  }
  throw new Error(`Server did not respond at ${url} within ${timeoutMs}ms`);
}

fs.mkdirSync(path.dirname(outFile), { recursive: true });

await run('npm', ['run', 'build'], { cwd: frontendRoot });

const preview = startPreview();
let previewExited = false;
preview.on('close', () => {
  previewExited = true;
});

try {
  await waitForHttpOk('http://127.0.0.1:4173/', 60000);
  await new Promise((r) => setTimeout(r, 1200));

  const browser = await chromium.launch();
  const page = await browser.newPage({
    viewport: { width: 1600, height: 1100 },
    deviceScaleFactor: 2,
  });
  await page.goto('http://127.0.0.1:4173/', { waitUntil: 'networkidle' });
  await page.waitForSelector('.react-flow__node', { timeout: 30000 });
  await new Promise((r) => setTimeout(r, 800));

  const panel = page.locator('[data-longevity-diagram="capture"]');
  await panel.screenshot({ path: outFile, type: 'png' });

  await browser.close();
  console.log(`Wrote ${outFile}`);
} finally {
  if (!previewExited) {
    preview.kill('SIGTERM');
  }
}
