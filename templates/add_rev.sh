#!/bin/bash
psdb new_rev ${default_toml_file} ${revname}
cat >"${revname}.toml" <<EOL
[worker]
nodes = "\"${nodes}\""
max_concurrent = ${max_concurrent}
run_on_file_host = true
run_on_file_host_pattern = '\\/net/(node\d{3})'
EOL
