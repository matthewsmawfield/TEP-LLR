#!/usr/bin/env node

const fs = require('fs');
const path = require('path');

function readJsonIfExists(filePath) {
    if (!fs.existsSync(filePath)) {
        return null;
    }
    const raw = fs.readFileSync(filePath, 'utf8');
    try {
        return JSON.parse(raw);
    } catch (error) {
        throw new Error(`Failed to parse JSON ${filePath}: ${error.message}`);
    }
}

function formatFixedNumber(value, decimals) {
    if (value === null || value === undefined) {
        return '';
    }
    const num = Number(value);
    if (!Number.isFinite(num)) {
        return '';
    }
    return num.toFixed(decimals);
}

function formatSignedNumber(value, decimals) {
    const fixed = formatFixedNumber(value, decimals);
    if (!fixed) {
        return '';
    }
    const num = Number(value);
    return `${num >= 0 ? '+' : ''}${fixed}`;
}

function formatCI(lower, upper, decimals) {
    const lo = formatSignedNumber(lower, decimals);
    const hi = formatSignedNumber(upper, decimals);
    if (!lo || !hi) {
        return '';
    }
    return `[${lo}, ${hi}]`;
}

function formatIntegerWithCommas(value) {
    if (value === null || value === undefined) {
        return '';
    }
    const num = Number(value);
    if (!Number.isFinite(num)) {
        return '';
    }
    return Math.round(num).toLocaleString('en-US');
}

function formatIntegerLatex(value) {
    const s = formatIntegerWithCommas(value);
    if (!s) {
        return '';
    }
    return s.replace(/,/g, '{,}');
}

function bayesFactorEvidenceLabel(bf) {
    if (bf === null || bf === undefined) {
        return '';
    }
    const num = Number(bf);
    if (!Number.isFinite(num) || num < 0) {
        return '';
    }
    if (num < 1) {
        return 'negative';
    }
    if (num < 3) {
        return 'barely worth mentioning';
    }
    if (num < 10) {
        return 'substantial';
    }
    if (num < 30) {
        return 'strong';
    }
    if (num < 100) {
        return 'very strong';
    }
    return 'decisive';
}

function formatPValueMantissaExp(pValue, mantissaDecimals = 1) {
    if (pValue === null || pValue === undefined) {
        return null;
    }
    const p = Number(pValue);
    if (!Number.isFinite(p) || p <= 0) {
        return null;
    }
    const exp = Math.floor(Math.log10(p));
    let mantissa = p / Math.pow(10, exp);
    let mantissaStr = mantissa.toFixed(mantissaDecimals);
    if (Number(mantissaStr) >= 10) {
        mantissa = mantissa / 10;
        mantissaStr = mantissa.toFixed(mantissaDecimals);
        return { mantissa: mantissaStr, exp: exp + 1 };
    }
    return { mantissa: mantissaStr, exp };
}

function formatPValueLatex(pValue, mantissaDecimals = 1) {
    const parts = formatPValueMantissaExp(pValue, mantissaDecimals);
    if (!parts) {
        return null;
    }
    return `${parts.mantissa} \\times 10^{${parts.exp}}`;
}

function formatScientificLatex(value, mantissaDecimals = 2) {
    if (value === null || value === undefined) {
        return null;
    }
    const num = Number(value);
    if (!Number.isFinite(num)) {
        return null;
    }
    if (num === 0) {
        return '0';
    }
    const sign = num < 0 ? '-' : '';
    const magnitude = Math.abs(num);
    let exp = Math.floor(Math.log10(magnitude));
    let mantissa = magnitude / Math.pow(10, exp);
    let mantissaStr = mantissa.toFixed(mantissaDecimals);
    if (Number(mantissaStr) >= 10) {
        mantissa = mantissa / 10;
        mantissaStr = mantissa.toFixed(mantissaDecimals);
        exp += 1;
    }
    return `${sign}${mantissaStr} \\times 10^{${exp}}`;
}

function formatPValueCell(pValue, options = {}) {
    const { includePrefix = false } = options;
    if (pValue === null || pValue === undefined) {
        return '';
    }
    const p = Number(pValue);
    if (!Number.isFinite(p)) {
        return '';
    }
    if (p >= 0.001) {
        const s = formatFixedNumber(p, 3);
        return includePrefix ? `p = ${s}` : s;
    }
    const latex = formatPValueLatex(p);
    if (!latex) {
        return '';
    }
    return includePrefix ? `p = ${latex}` : latex;
}

function formatPValueHtmlCell(pValue) {
    if (pValue === null || pValue === undefined) {
        return '';
    }
    const p = Number(pValue);
    if (!Number.isFinite(p) || p <= 0) {
        return '';
    }
    if (p >= 0.001) {
        return formatFixedNumber(p, 3);
    }
    const latex = formatPValueLatex(p);
    if (!latex) {
        return '';
    }
    return `$${latex}$`;
}

function safeGet(obj, pathParts) {
    let cur = obj;
    for (const part of pathParts) {
        if (cur === null || cur === undefined) {
            return undefined;
        }
        cur = cur[part];
    }
    return cur;
}

function createInjectionContext() {
    const outputsDir = path.join(__dirname, '..', 'results', 'outputs');

    // Load TEP-LLR pipeline outputs used by manuscript placeholder injection.
    const step003 = readJsonIfExists(path.join(outputsDir, 'step_003_statistical_analysis.json'));
    const step004 = readJsonIfExists(path.join(outputsDir, 'step_004_detection_analysis_advanced.json'));
    const step006 = readJsonIfExists(path.join(outputsDir, 'step_006_multi_ephemeris_comparison.json'));
    const step016 = readJsonIfExists(path.join(outputsDir, 'step_016_bayesian_analysis.json'));
    const step017 = readJsonIfExists(path.join(outputsDir, 'step_017_leverage_diagnostics.json'));
    const step029 = readJsonIfExists(path.join(outputsDir, 'step_029_station_power_analysis.json'));
    const step030 = readJsonIfExists(path.join(outputsDir, 'step_030_hardware_epoch_analysis.json'));

    const ctx = {
        llr: {
            step002: {},
            step016: {},
            step031: {},
            step032: {},
            step018: {},
            step005: {},
            tep: {}
        }
    };

    // Step 002: Primary statistical analysis
    if (step003) {
        const results = step003.analysis_results || step003.regression_metrics || {};
        ctx.llr.step002 = {
            n_observations: formatIntegerWithCommas(results.n_obs || step003.outlier_cleaning?.n_cleaned),
            eta_ols_sci: results.eta ?? step003.eta_ols,
            eta_err_sci: results.eta_error ?? step003.eta_ols_error,
            snr: formatFixedNumber(step003.snr, 2),
            sigma: formatFixedNumber(step003.snr, 2)
        };
    } else {
        console.warn('⚠️  Missing step_003_statistical_analysis.json; primary detection values not injected.');
    }

    // Step 016: Bayesian analysis
    if (step016) {
        const bayes = step016.bayesian_summary || {};
        const bf = bayes.bayes_factor_savage_dickey;
        ctx.llr.step016 = {
            bayes_factor_llr_sci: bf || 0,
            bayes_factor_llr_val: formatScientificLatex(bf),
            evidence_strength: bayesFactorEvidenceLabel(bf)
        };
    }

    // Step 032: Hardware epoch analysis
    if (step030) {
        const epochs = [
            ...(step030.epoch_fits?.grasse_epochs || []),
            ...(step030.epoch_fits?.apo_epochs || []),
        ];
        const modern = epochs.find(e => e.epoch_name === 'Grasse-III') || {};
        ctx.llr.step032 = {
            modern_snr: formatFixedNumber(modern.snr, 2),
            modern_eta: modern.eta
        };
    }

    // Step 031: Station power analysis
    if (step029) {
        const stations = step029.per_station_power?.stations || [];
        const apo = stations.find(s => s.station === 'APO') || {};
        const grasse = stations.find(s => s.station === 'Grasse') || {};
        const pw = step029.precision_weighted_regression || {};
        
        ctx.llr.step031 = {
            apo_eta_sci: apo.eta_obs,
            apo_err_sci: apo.eta_err_obs,
            apo_snr: formatFixedNumber(apo.snr_observed, 1),
            apo_n: formatIntegerWithCommas(apo.n_obs),
            grasse_eta_sci: grasse.eta_obs,
            grasse_snr: formatFixedNumber(grasse.snr_observed, 1),
            pw_eta_sci: pw.eta_precision_weighted,
            pw_snr: formatFixedNumber(pw.snr, 2),
            cross_station_r: formatFixedNumber(step029.cross_station_validation?.prediction_r, 4)
        };
    } else {
        console.warn('⚠️  Missing step_029_station_power_analysis.json; station values not injected.');
    }

    // Step 018: Leverage diagnostics
    if (step017) {
        const cook = step017.conclusion?.formal_cooks_d_excision || {};
        const summary = step017.summary || {};
        
        ctx.llr.step018 = {
            theilsen_eta_sci: summary.full_sample_eta_theilsen,
            cook_eta_sci: cook.eta_clean_ols,
            cook_err_sci: cook.eta_clean_se,
            cook_snr: formatFixedNumber(cook.eta_clean_snr, 2),
            n_high_leverage: formatIntegerWithCommas(step017.leverage_statistics?.n_high_leverage)
        };
    } else {
        console.warn('⚠️  Missing step_017_leverage_diagnostics.json; leverage values not injected.');
    }

    // Step 005: DE430 cross-validation
    if (step006) {
        const de430 = step006.comparisons?.DE430 || {};
        ctx.llr.step005 = {
            de430_eta_sci: de430.eta,
            de430_err_sci: de430.eta_error,
            de430_snr: formatFixedNumber(de430.snr, 2)
        };
    } else {
        console.warn('⚠️  Missing step_006_multi_ephemeris_comparison.json; DE430 values not injected.');
    }

    // Backward compatibility: map to old tep.* paths
    ctx.tep = {
        table12: {},
        table12c: {},
        table12d: {},
        meta: {},
        step114: {}
    };

    return ctx;
}

function injectPlaceholders(template, context) {
    const unresolved = new Set();
    const replaced = template.replace(/\{\{\s*([^}]+?)\s*\}\}/g, (match, expr) => {
        const pathParts = expr.split('.').map(s => s.trim()).filter(Boolean);
        const value = safeGet(context, pathParts);
        if (value === undefined || value === null) {
            unresolved.add(expr);
            return match;
        }
        return String(value);
    });
    if (unresolved.size > 0) {
        throw new Error(`Unresolved placeholders (${unresolved.size}): ${Array.from(unresolved).join(', ')}`);
    }
    return replaced;
}

async function buildStaticSite() {
    console.log('🔨 Building static site...');
    
    try {
        // Clean dist directory
        const distDir = path.join(__dirname, 'dist');
        if (fs.existsSync(distDir)) {
            console.log('🧹 Cleaning dist directory...');
            fs.rmSync(distDir, { recursive: true, force: true });
        }
        
        // Read the manifest
        const manifestPath = path.join(__dirname, 'manifest.json');
        const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
        
        // Read the base index.html
        const indexPath = path.join(__dirname, 'index.html');
        let indexContent = fs.readFileSync(indexPath, 'utf8');
        
        const injectionContext = createInjectionContext();

        // Build the component content
        let componentsHtml = '';
        
        // Sort sections by order and load each component
        const sortedSections = manifest.sections.sort((a, b) => a.order - b.order);
        
        for (const section of sortedSections) {
            console.log(`📄 Loading section: ${section.title}`);
            
            const componentPath = path.join(__dirname, 'components', section.file);
            
            if (fs.existsSync(componentPath)) {
                const componentHtmlRaw = fs.readFileSync(componentPath, 'utf8');
                const componentHtml = injectPlaceholders(componentHtmlRaw, injectionContext);
                
                // Wrap component in section container (matching the dynamic loader)
                componentsHtml += `
                <section id="${section.id}" class="manuscript-section" data-section="${section.title}">
                    ${componentHtml}
                </section>`;
            } else {
                throw new Error(`Component not found: ${componentPath}`);
            }
        }
        
        // Replace the dynamic content with static content
        const staticContent = indexContent
            .replace(
                /<div id="loading".*?<\/div>\s*<div id="manuscript-content".*?<\/div>/s,
                () => `<div id="manuscript-content">${componentsHtml}</div>`
            )
            .replace(/<div id="loading"[^>]*>[\s\S]*?<\/div>/g, '')
            .replace(/<div id="manuscript-content"[^>]*style="display:\s*none;"[^>]*>[\s\S]*?<\/div>/g, '')
            .replace(
                /<!-- Component Loading Script -->[\s\S]*?<\/script>/g,
                `<!-- Static build - components pre-loaded -->\n    <script>\n        document.addEventListener("DOMContentLoaded", function() {\n            if (window.MathJax && window.MathJax.typesetPromise) {\n                window.MathJax.typesetPromise().catch(function (err) {\n                    console.error("MathJax error:", err.message);\n                });\n            }\n        });\n    </script>`
            )
            .replace(
                '<main id="main-content" role="main">',
                '<!-- This is a statically built version for SEO/deployment -->\n    <main id="main-content" role="main">'
            );
        
        // Create dist directory if it doesn't exist
        if (!fs.existsSync(distDir)) {
            fs.mkdirSync(distDir, { recursive: true });
        }
        
        // Write the built file
        const outputPath = path.join(distDir, 'index.html');
        fs.writeFileSync(outputPath, staticContent, 'utf8');
        
        // Copy necessary static assets to dist
        // EXPLICITLY COPY STYLES DIRECTORY
        const assetDirs = ['public', 'figures', 'data', 'styles'];
        for (const assetDir of assetDirs) {
            const srcPath = path.join(__dirname, assetDir);
            const destPath = path.join(distDir, assetDir);
            
            if (fs.existsSync(srcPath)) {
                console.log(`📁 Copying ${assetDir}/`);
                copyRecursiveSync(srcPath, destPath);
            }
        }

        // Copy simple-styles.css if it exists
        const simpleStylesSrc = path.join(__dirname, 'simple-styles.css');
        const simpleStylesDest = path.join(distDir, 'simple-styles.css');
        if (fs.existsSync(simpleStylesSrc)) {
            fs.copyFileSync(simpleStylesSrc, simpleStylesDest);
            console.log('📁 Copied simple-styles.css');
        }
        
        // Copy figures from results/figures/ to dist/figures/ (main figure source)
        const resultsFiguresPath = path.join(__dirname, '..', 'results', 'figures');
        const distFiguresPath = path.join(distDir, 'figures');
        const distPublicFiguresPath = path.join(distDir, 'public', 'figures');
        if (fs.existsSync(resultsFiguresPath)) {
            console.log('📁 Copying results/figures/ → dist/figures/');
            copyRecursiveSync(resultsFiguresPath, distFiguresPath);
            console.log('📁 Copying results/figures/ → dist/public/figures/');
            copyRecursiveSync(resultsFiguresPath, distPublicFiguresPath);
        }
        
        // Copy manifest.json for reference
        fs.copyFileSync(manifestPath, path.join(distDir, 'manifest.json'));
        
        // Copy .nojekyll to dist root for GitHub Pages
        const nojekyllSrc = path.join(__dirname, 'public', '.nojekyll');
        const nojekyllDest = path.join(distDir, '.nojekyll');
        if (fs.existsSync(nojekyllSrc)) {
            fs.copyFileSync(nojekyllSrc, nojekyllDest);
            console.log('📁 Copied .nojekyll to dist root');
        } else {
            // Create it if it doesn't exist
            fs.writeFileSync(nojekyllDest, '');
            console.log('📁 Created .nojekyll in dist root');
        }

        // Copy robots.txt and sitemap.xml to dist root
        const rootFiles = ['404.html', 'robots.txt', 'sitemap.xml', 'CNAME', '29c6507763d2303d801cc8ed89d39f88.txt', 'favicon.ico'];
        for (const file of rootFiles) {
            const src = path.join(__dirname, 'public', file);
            const dest = path.join(distDir, file);
            if (fs.existsSync(src)) {
                fs.copyFileSync(src, dest);
                console.log(`📁 Copied ${file} to dist root`);
            }
        }
        
        // Copy citation files to dist root
        const citationFiles = ['CITATION.cff', 'CITATION.bib', 'citation.json', 'codemeta.json', 'README.md'];
        const optionalCitationFiles = new Set(['LICENSE']);
        const projectRoot = path.join(__dirname, '..');
        
        for (const file of [...citationFiles, ...optionalCitationFiles]) {
            const candidates = [
                path.join(projectRoot, file),
                path.join(__dirname, file),
            ];
            const src = candidates.find(candidate => fs.existsSync(candidate));
            const dest = path.join(distDir, file);
            if (src) {
                fs.copyFileSync(src, dest);
                console.log(`📁 Copied ${file} to dist root`);
            } else if (!optionalCitationFiles.has(file)) {
                console.warn(`⚠️  Missing citation file: ${file}`);
            }
        }
        
        // Copy .well-known directory
        const wellKnownSrc = path.join(__dirname, 'public', '.well-known');
        const wellKnownDest = path.join(distDir, '.well-known');
        if (fs.existsSync(wellKnownSrc)) {
            if (!fs.existsSync(wellKnownDest)) {
                fs.mkdirSync(wellKnownDest, { recursive: true });
            }
            const files = fs.readdirSync(wellKnownSrc);
            for (const file of files) {
                fs.copyFileSync(path.join(wellKnownSrc, file), path.join(wellKnownDest, file));
            }
            console.log(`📁 Copied .well-known/ to dist`);
        }
        
        // Generate markdown version
        console.log('📝 Generating markdown version...');
        const { HTMLToMarkdownConverter } = require('./html-to-markdown.js');
        const converter = new HTMLToMarkdownConverter();
        await converter.convertSiteToMarkdown();
        
        console.log('✅ Static site built successfully!');
        console.log(`📁 Output: ${outputPath}`);
        console.log(`📊 Generated ${manifest.sections.length} sections`);
        console.log('🚀 Ready for deployment');
        
    } catch (error) {
        console.error('❌ Build failed:', error);
        process.exit(1);
    }
}

// Helper function to copy directories recursively
function copyRecursiveSync(src, dest) {
    const exists = fs.existsSync(src);
    const stats = exists && fs.statSync(src);
    const isDirectory = exists && stats.isDirectory();
    
    if (isDirectory) {
        if (!fs.existsSync(dest)) {
            fs.mkdirSync(dest, { recursive: true });
        }
        fs.readdirSync(src).forEach(childItemName => {
            copyRecursiveSync(
                path.join(src, childItemName),
                path.join(dest, childItemName)
            );
        });
    } else {
        // Filter out raw data files
        if (src.endsWith('.csv') || src.endsWith('.dat') || src.endsWith('.nc')) {
            return;
        }
        fs.copyFileSync(src, dest);
    }
}

// Run if called directly
if (require.main === module) {
    buildStaticSite();
}

module.exports = { buildStaticSite };
