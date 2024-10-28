import numpy as np

str1 = "1"
str2 = "2"
str3 = "3"
str4 = "4"
guess1 = 121
guess2 = 57
guess3 = 31
guess4 = 205

strs = [str1, str2, str3, str4]
guesses = [guess1, guess2, guess3, guess4]

print(strs[np.argmin(guesses)])