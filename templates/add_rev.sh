#!/bin/bash

mkdir -p "${ps_dir}"
cd "${ps_dir}"
cp "${projectDir}/configs/pspipe_toml_templates/default.toml" .
cp "${projectDir}/configs/pspipe_toml_templates/eor_bins_hba.parset" .
cp "${projectDir}/configs/pspipe_toml_templates/ps_config_hba.parset" .
cp "${projectDir}/configs/pspipe_toml_templates/flagger.parset" .
cp "${projectDir}/configs/pspipe_toml_templates/gpr_config_hba.parset" .
cp "${projectDir}/configs/pspipe_toml_templates/gpr_config_v.parset" .
cp "${projectDir}/configs/pspipe_toml_templates/gpr_ml_config_2023.parset" .
cp "${projectDir}/configs/pspipe_toml_templates/flagger_pre_combine.parset" .

psdb new_rev ${ps_dir}/default.toml ${revname}
cat >"${revname}.toml" <<EOL
default_settings = "${ps_dir}/default.toml"
data_dir = "${ps_dir}"
[worker]
nodes = \"${nodes}\"
max_concurrent = ${max_concurrent}
run_on_file_host = true
run_on_file_host_pattern = '\\/net/(node\\d{3})'
env_file='/home/users/mertens/.activate_pspipe_dev.sh'
[merge_ms]
data_col = "${merge_data_column}"
#aoflagger_strategy = '/home/users/mertens/projects/NCP/nights_np5_red1/lofar-sens2.lua'
[image]
data_col = "${image_data_col}"
channels_out = 'every3'
name="${revname}"
[power_spectra]
eor_bin_list = "${ps_dir}/eor_bins_hba.parset"
ps_config = "${ps_dir}/ps_config_hba.parset"
flagger = "${ps_dir}/flagger.parset"
[gpr]
config_i = "${ps_dir}/gpr_config_hba.parset"
config_v = "${ps_dir}/gpr_config_v.parset"
[ml_gpr]
name = 'eor_vae_2023'
config = "${ps_dir}/gpr_ml_config_2023.parset"
[combine]
pre_flag = "${ps_dir}/flagger_pre_combine.parset"
EOL
