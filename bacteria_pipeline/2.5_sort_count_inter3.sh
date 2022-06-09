#!/bin/bash

echo 'running... ' $(basename $0)


DIR=$1

sort ${DIR}inter3_replaced | uniq -c > ${DIR}inter4_replaced_sort_count


