#!/bin/sh
# Regenerate static site and Lucknow markdown without invoking npm (avoids zsh
# compdef / completion noise when npm is hooked into the interactive shell).
set -eu
ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
cd "$ROOT/site"
node build.js
node html-to-markdown.js
