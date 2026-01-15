# R Path Configuration

## Issue
Cursor/VS Code cannot find R for help, package installation, etc.

## Solution

### Step 1: Install R (if not already installed)

**Option A: Install via Homebrew (Recommended)**
```bash
brew install r
```

**Option B: Download from CRAN**
1. Visit: https://cran.r-project.org/bin/macosx/
2. Download and install the appropriate version for your Mac

### Step 2: Find Your R Installation Path

Run the helper script:
```bash
bash find_r.sh
```

Or manually check these common locations:
- `/Library/Frameworks/R.framework/Resources/bin/R` (CRAN installer)
- `/usr/local/bin/R` (Homebrew on Intel Mac)
- `/opt/homebrew/bin/R` (Homebrew on Apple Silicon)
- `/usr/bin/R` (system R, if available)

### Step 3: Update Settings

Edit `.vscode/settings.json` and set the `r.rpath.mac` value to your R path:

```json
{
  "r.rpath.mac": "/opt/homebrew/bin/R"
}
```

Replace `/opt/homebrew/bin/R` with the actual path found in Step 2.

### Step 4: Reload Cursor/VS Code

After updating the settings, reload the window (Cmd+Shift+P → "Reload Window") or restart Cursor/VS Code.
