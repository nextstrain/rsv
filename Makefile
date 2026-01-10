# Makefile for Nextstrain RSV pipeline

# Generate YAML config from CUE
config/configfile.yaml: config/configfile.cue
	@echo "Generating YAML from CUE..."
	micromamba run -n cue cue export --force $< --out yaml --outfile $@
	@echo "✓ Generated $@"

# Validate CUE syntax
.PHONY: validate-config
validate-config:
	@echo "Validating CUE syntax..."
	micromamba run -n cue cue vet config/configfile.cue
	@echo "✓ CUE validation passed"

# Format CUE file
.PHONY: format-config
format-config:
	@echo "Formatting CUE file..."
	micromamba run -n cue cue fmt config/configfile.cue
	@echo "✓ CUE file formatted"

# Test: validate and export
.PHONY: test-config
test-config: validate-config config/configfile.yaml
	@echo "✓ Config test passed"

.PHONY: help
help:
	@echo "Nextstrain RSV Pipeline - Available targets:"
	@echo ""
	@echo "Config management:"
	@echo "  config/configfile.yaml    Generate YAML from CUE source"
	@echo "  validate-config           Validate CUE syntax"
	@echo "  format-config             Format CUE file"
	@echo "  test-config               Validate and export CUE"
	@echo ""
	@echo "Note: Edit config/configfile.cue (not .yaml) to modify config"
