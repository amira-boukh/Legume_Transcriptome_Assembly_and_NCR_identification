#!/usr/bin/env python3
"""Print True if ranking_debug.json has mean pLDDT equal or higher tan 70, else False"""
import json
import sys

if len(sys.argv) != 2:
    print("False")
    sys.exit(0)

with open(sys.argv[1], "r", encoding="utf-8") as handle:
    data = json.load(handle)

scores = data.get("plddts") or data.get("plddt") or data.get("ranking_confidences") or {}
if isinstance(scores, dict):
    value = max(scores.values()) if scores else 0.0
elif isinstance(scores, list):
    value = max(scores) if scores else 0.0
elif isinstance(scores, (int, float)):
    value = float(scores)
else:
    value = 0.0

print("True" if value >= 70 else "False")
