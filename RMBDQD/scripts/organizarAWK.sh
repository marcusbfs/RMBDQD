#!/bin/bash

cat * | awk  ' BEGIN{printf "%s\t%s\n", "rho","E"} 
/rho|E\s+/ {
if($1=="E"){printf "\t%s\n", $3} 
else if ($1=="rho")
{printf "%.2f",$3}}' > "$1"

sed -i '/^0.*\.0\+$/d' "$1"
