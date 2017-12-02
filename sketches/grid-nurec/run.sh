mkdir -p par data cool img
parallel bash job.sh {/} ::: par/*
