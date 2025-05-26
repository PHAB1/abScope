awk 'BEGIN{FS="\t"; OFS="\t"} {
  dupcount_value = 1
  if (match($0, /DUPCOUNT=([0-9]+)/, arr)) {
    dupcount_value = arr[1]
  }
  print $0, dupcount_value
}' all_phage_airr.tsv > all_treat_airr_dup.tsv
