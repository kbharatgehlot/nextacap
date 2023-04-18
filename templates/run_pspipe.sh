#!/bin/bash

cd "${ps_dir}"

psdb add_obs ${toml_file} ${obsid} -m ${msfiles}
obs="${obsid}"

if ${merge_ms}; then
    echo "Merging Ms files"
    if ${delay_flag}; then
        echo "Using delay flagger"
        pspipe merge_ms,delay_flagger ${toml_file} ${obsid}
    else
        pspipe merge_ms ${toml_file} ${obsid}
    fi

    obs="${obsid}_flagged"
fi
echo "making image cube"
pspipe image,gen_vis_cube ${toml_file} \${obs}

if ${ml_gpr}; then
    echo "running ML_GPR"
    pspipe run_ml_gpr ${toml_file} \${obs}
else
    echo "running GPR"
    pspipe run_gpr ${toml_file} \${obs}
fi
