primeiro remapeamos
python3 ../../scripts/script_analysis/map_ids2orig.py /home/glicina/pbarros/abScope/results/quality_trimmed/nano_err_treat/qc.meta nano_treat_error.tsv

Rename '*' e '|'
bash ../../scripts/script_analysis/recursive_id_rename.sh

define clones column
bash ../../scripts/script_analysis/recursive_define_clones.sh . define_clones/

Join all
python3 ../../scripts/script_analysis/recursiveAnnotate.py -i define_clones/ -o .

Insere coluna duplicados com base no nome da sequencia DUPCOUNT=X ; substitui head para ser "duplicate_count"
bash ../../scripts/script_analysis/dupcount_col.sh
