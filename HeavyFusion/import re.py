import re

s = "100e-2 2143e+10 1234e-3 123e4"
pattern = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"

matches = re.findall(pattern, s)
print(matches)
# 转为浮点数
matches = [float(match) for match in matches]
print(matches)