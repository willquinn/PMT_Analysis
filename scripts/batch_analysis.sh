# USER INPUTS:
# ============
while getopts :t: option;
do
  case "${option}" in
      t )
        tag=${OPTARG}
        ;;
  esac
done

root_dir="/unix/nemo4/PMT_He_Study_nemo4"
script_log_dir=batch_log_$tag
mkdir $root_dir/batch_files/$script_log_dir
# mkdir $root_dir/batch_files/$script_log_dir/scripts
mkdir $root_dir/batch_files/$script_log_dir/logs

# Files found  automatically
ls $root_dir/data/ROOT_files/*.root > $root_dir/batch_files/$script_log_dir/filenames.ascii
data_file_list=`cat $root_dir/batch_files/$script_log_dir/filenames.ascii`
#
config_file="/home/wquinn/pmt_analysis/PMT_Analysis/config_files/pmt_permeation_config_file.txt"
source_code="/home/wquinn/pmt_analysis/PMT_Analysis/scripts/root_reader.py"

ifile=0
for file in $data_file_list;
do
  executable_plus_arguments="python3 $source_code -i $file -c $config_file -o /home/wquinn/pmt_analysis/output.txt"

  interactive_command="$executable_plus_arguments >& $root_dir/batch_files/$script_log_dir/logs/batch_file_$ifile.log"
  echo $interactive_command
  $interactive_command
done