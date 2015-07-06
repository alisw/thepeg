SUBDIRS= ThePEG Pythia7 TheP8I test Ariadne

check-local:
	@find . -maxdepth 2 -type f -name '*.log' -print0 | xargs -0 -e grep -n -e 'The following exception' || echo "No exceptions found in log files."

EXTRA_DIST = GUIDELINES
