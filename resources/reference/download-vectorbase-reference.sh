#!/bin/bash
set -euo pipefail

# Default values
GENOME_NAME="AgambiaePEST"
VERSION="66"
STRIP_STRING="AgamP4_"
OUTPUT_DIR="./resources/reference"

# vectorbase.org/common/downloads/... was retired (returns a 200/404 HTML
# "Page Not Found" page for every release, on vectorbase.org and other
# VEuPathDB sites alike) when VectorBase moved to a JS-rendered downloads
# app. The same files are still served, at the same release-numbered path,
# from this legacy host.
DOWNLOAD_HOST="https://legacy.vectorbase.org"

# Function to display usage information
usage() {
    echo "Usage: $0 [GENOME_NAME] [VERSION] [STRIP_STRING]"
    echo "  GENOME_NAME   : Name of the genome (default: AgambiaePEST)"
    echo "  VERSION       : Version number (default: 66)"
    echo "  STRIP_STRING  : String to strip from sequence headers (default: AgamP4_)"
    echo ""
    echo "Example: $0 AgambiaePEST 66 AgamP4_"
    exit 1
}

# Parse arguments
[ $# -ge 1 ] && GENOME_NAME=$1
[ $# -ge 2 ] && VERSION=$2
[ $# -ge 3 ] && STRIP_STRING=$3

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Extract the short name from the STRIP_STRING (removing trailing underscore if present)
SHORT_NAME=${STRIP_STRING%_}
[ -z "$SHORT_NAME" ] && SHORT_NAME=$GENOME_NAME

echo "Downloading genome and annotation files for $GENOME_NAME (version $VERSION)"
echo "Will strip '$STRIP_STRING' from sequence headers"

# Downloads a URL and sanity-checks the result before trusting it: --fail
# makes curl itself error out on a 4xx/5xx response instead of quietly
# emitting the server's HTML error page as if it were the file, and the
# leading '>' / '##gff-version' check below catches the case where a 200
# response still isn't the expected format (e.g. an HTML app shell).
fetch() {
    local url=$1 out=$2 expect_prefix=$3 status

    # Temporarily disable -e so a failed curl doesn't abort the script before
    # we can report which URL failed and clean up the partial output file.
    set +e
    curl -fL --progress-bar "$url" | sed "s/$STRIP_STRING//g" > "$out"
    status=${PIPESTATUS[0]}
    set -e

    if [ "$status" -ne 0 ]; then
        echo "ERROR: download failed (curl exit $status) for $url" >&2
        rm -f "$out"
        exit 1
    fi
    if ! head -c "${#expect_prefix}" "$out" | grep -qF "$expect_prefix"; then
        echo "ERROR: $out does not look like the expected file (doesn't start with '$expect_prefix')." >&2
        echo "       $url likely did not return real data - check the URL/host/release version." >&2
        rm -f "$out"
        exit 1
    fi
}

# Download and process the FASTA file
echo "Downloading FASTA file..."
FASTA_URL="$DOWNLOAD_HOST/common/downloads/release-$VERSION/$GENOME_NAME/fasta/data/VectorBase-${VERSION}_${GENOME_NAME}_Genome.fasta"
fetch "$FASTA_URL" "$OUTPUT_DIR/$SHORT_NAME.fa" ">"

# Download and process the GFF file
echo "Downloading GFF file..."
GFF_URL="$DOWNLOAD_HOST/common/downloads/release-$VERSION/$GENOME_NAME/gff/data/VectorBase-${VERSION}_${GENOME_NAME}.gff"
fetch "$GFF_URL" "$OUTPUT_DIR/$SHORT_NAME.gff" "##gff-version"

echo "Download complete. Files saved to $OUTPUT_DIR/"
echo "  - $OUTPUT_DIR/$SHORT_NAME.fa"
echo "  - $OUTPUT_DIR/$SHORT_NAME.gff"
