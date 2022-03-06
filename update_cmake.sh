#!/bin/bash


getopts ":p:" opt
  case $opt in
    p)
      ## echo "-p was triggered, Parameter: $OPTARG" >&2
      for d in ./*/; do
	      awk '/Include universal/ { print; print "include_directories(\"'$OPTARG'\")"; next }1' $d/CMakeLists.txt > tmp_CMakeLists.txt
	      # Need a temp file
	      cat tmp_CMakeLists.txt > $d/CMakeLists.txt
	      rm tmp_CMakeLists.txt
      done
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      echo "Use -p <universal directory>" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
