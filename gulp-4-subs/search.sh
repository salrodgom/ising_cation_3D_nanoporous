#!/bin/bash
grep " $1 " OUTSOD* | grep " $2 " | grep " $3 " | grep " $4"
