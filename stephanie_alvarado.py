from scipy import stats
from scipy.stats import t as t_dist
from scipy.stats import chi2

from abtesting_test import *

# You can comment out these lines! They are just here to help follow along to the tutorial.
# print(t_dist.cdf(-2, 20)) # should print .02963
# print(t_dist.cdf(2, 20)) # positive t-score (bad), should print .97036 (= 1 - .2963)
#
# print(chi2.cdf(23.6, 12)) # prints 0.976
# print(1 - chi2.cdf(23.6, 12)) # prints 1 - 0.976 = 0.023 (yay!)

# TODO: Fill in the following functions! Be sure to delete "pass" when you want to use/run a function!
# NOTE: You should not be using any outside libraries or functions other than the simple operators (+, **, etc)
# and the specifically mentioned functions (i.e. round, cdf functions...)

def slice_2D(list_2D, start_row, end_row, start_col, end_col):
    '''
    Splices a the 2D list via start_row:end_row and start_col:end_col
    :param list: list of list of numbers
    :param nums: start_row, end_row, start_col, end_col
    :return: the spliced 2D list (ending indices are exclsive)
    '''
    to_append = []
    for l in range(start_row, end_row):
        to_append.append(list_2D[l][start_col:end_col])

    return to_append

def get_avg(nums):
    '''
    Helper function for calculating the average of a sample.
    :param nums: list of numbers
    :return: average of list
    '''
    total = 0
    for i in nums:
        total+=i
    return total/len(nums)


def get_stdev(nums):
    '''
    Helper function for calculating the standard deviation of a sample.
    :param nums: list of numbers
    :return: standard deviation of list
    '''
    average = get_avg(nums)
    sum = 0
    for i in nums:
        sum+=((i-average)**(2))
    return (sum/(len(nums)-1))**(1/2)

def get_standard_error(a, b):
    '''
    Helper function for calculating the standard error, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: standard error of a and b (see studio 6 guide for this equation!)
    '''
    stda_squared = get_stdev(a)**2
    stdb_squared = get_stdev(b)**2
    sum = (stda_squared/len(a)) + (stdb_squared/len(b))
    return sum**(1/2)

def get_2_sample_df(a, b):
    '''
    Calculates the combined degrees of freedom between two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: integer representing the degrees of freedom between a and b (see studio 6 guide for this equation!)
    HINT: you can use Math.round() to help you round!
    '''
    standard_error_4 = get_standard_error(a, b)**(4)
    len_a = len(a)
    len_b = len(b)
    part_1 = (((get_stdev(a)**2)/len_a)**2)/(len_a-1)
    part_2 = (((get_stdev(b)**2)/len_b)**2)/(len_b-1)
    denominator = part_1 + part_2
    return round(standard_error_4/denominator)

def get_t_score(a, b):
    '''
    Calculates the t-score, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: number representing the t-score given lists a and b (see studio 6 guide for this equation!)
    '''
    t_score = (get_avg(a) - get_avg(b)) / get_standard_error(a, b)
    if t_score > 0:
        t_score = (-1)*t_score

    return t_score

def perform_2_sample_t_test(a, b):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates a p-value by performing a 2-sample t-test, given two lists of numbers.
    :param a: list of numbers
    :param b: list of numbers
    :return: calculated p-value
    HINT: the t_dist.cdf() function might come in handy!
    '''
    return t_dist.cdf(get_t_score(a, b), get_2_sample_df(a,b))



# [OPTIONAL] Some helper functions that might be helpful in get_expected_grid().

def row_sum(observed_grid, ele_row):
    '''
    observed_grid: 2d list (row, col)
    ele_row: row number
    '''
    sum = 0
    for i in observed_grid[ele_row]:
        sum+=i
    return sum

def col_sum(observed_grid, ele_col):
    '''
    observed_grid: 2d list (row,col)
    ele_col: col number
    '''
    sum = 0
    #for every row at a specific column number
    for i in range(len(observed_grid)):
        sum+=observed_grid[i][ele_col]
    return sum

def total_sum(observed_grid):
    '''
        observed_grid: 2d list (row, col)
    '''
    total_sum = 0
    for row in range(len(observed_grid)):
        summed_row = row_sum(observed_grid, row)
        total_sum+=summed_row
    return total_sum

def calculate_expected(row_sum, col_sum, tot_sum):
    return (row_sum*col_sum)/tot_sum

def get_expected_grid(observed_grid):
    '''
    Calculates the expected counts, given the observed counts.
    ** DO NOT modify the parameter, observed_grid. **
    :param observed_grid: 2D list of observed counts
    :return: 2D list of expected counts
    HINT: To clean up this calculation, consider filling in the optional helper functions below!
    '''
    rows = len(observed_grid)
    cols = len(observed_grid[0])
    expected_grid = [[0]*cols for _ in range(rows)]
    total = total_sum(observed_grid)
    for row in range(rows):
        summed_row = row_sum(observed_grid, row)
        for col in range(cols):
            summed_col = col_sum(observed_grid, col)
            expected_value = calculate_expected(summed_row, summed_col, total)
            expected_grid[row][col] = expected_value
    return expected_grid


def df_chi2(observed_grid):
    '''
    Calculates the degrees of freedom of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: degrees of freedom of expected counts (see studio 6 guide for this equation!)
    '''
    rows = len(observed_grid)-1
    cols = len(observed_grid[0])-1
    return rows*cols

def chi2_value(observed_grid):
    '''
    Calculates the chi^2 value of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: associated chi^2 value of expected counts (see studio 6 guide for this equation!)
    '''
    chi2 = 0
    expected_grid = get_expected_grid(observed_grid)
    rows = len(expected_grid)
    cols = len(expected_grid[0])
    for row in range(rows):
        for col in range(cols):
            numerator = (observed_grid[row][col] - expected_grid[row][col])**2
            chi2+= (numerator/expected_grid[row][col])
    return chi2

def perform_chi2_homogeneity_test(observed_grid):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates the p-value by performing a chi^2 test, given a list of observed counts
    :param observed_grid: 2D list of observed counts
    :return: calculated p-value
    HINT: the chi2.cdf() function might come in handy!
    '''
    #TODO: fill me in!
    return 1-chi2.cdf(chi2_value(observed_grid), df_chi2(observed_grid))

# These commented out lines are for testing your main functions.
# Please uncomment them when finished with your implementation and confirm you get the same values :)
def data_to_num_list(s):
  '''
    Takes a copy and pasted row/col from a spreadsheet and produces a usable list of nums.
    This will be useful when you need to run your tests on your cleaned log data!
    :param str: string holding data
    :return: the spliced list of numbers
    '''
  return list(map(float, s.split()))


a_times ="""6902
315872
5687
14061
56734
14728
106835
102587
45000
127989
51378
15091"""

b_times ="""17604
10861
8669
4227
10944
94854
65311
29781
20667
45918
21289"""

a_return = """10 2"""

b_return = """8	3"""

complete_time_a = data_to_num_list(a_times)
complete_time_b = data_to_num_list(b_times)

print('TSCORE:' + str(get_t_score(complete_time_a, complete_time_b)))
print('P-VALUE: '+ str(perform_2_sample_t_test(complete_time_a, complete_time_b)))

a_returns = data_to_num_list(a_return)
b_returns = data_to_num_list(b_return)
returns_grid = [a_returns, b_returns]

print('chi_2 return' + str(chi2_value(returns_grid))) # this should be 4.103536
print('p-value return: ' + str(perform_chi2_homogeneity_test(returns_grid))) # this should be .0427939
#
# # t_test 1:
# a_t1_list = data_to_num_list(a1)
# b_t1_list = data_to_num_list(b1)
# print(get_t_score(a_t1_list, b_t1_list)) # this should be -129.500
# print(perform_2_sample_t_test(a_t1_list, b_t1_list)) # this should be 0.0000
# # why do you think this is? Take a peek at a1 and b1 in abtesting_test.py :)
#
# # t_test 2:
# a_t2_list = data_to_num_list(a2)
# b_t2_list = data_to_num_list(b2)
# print(get_t_score(a_t2_list, b_t2_list)) # this should be -1.48834
# print(perform_2_sample_t_test(a_t2_list, b_t2_list)) # this should be .082379
#
# # t_test 3:
# a_t3_list = data_to_num_list(a3)
# b_t3_list = data_to_num_list(b3)
# print(get_t_score(a_t3_list, b_t3_list)) # this should be -2.88969
# print(perform_2_sample_t_test(a_t3_list, b_t3_list)) # this should be .005091
#
# # chi2_test 1:
# a_c1_list = data_to_num_list(a_count_1)
# b_c1_list = data_to_num_list(b_count_1)
# c1_observed_grid = [a_c1_list, b_c1_list]
# print(chi2_value(c1_observed_grid)) # this should be 4.103536
# print(perform_chi2_homogeneity_test(c1_observed_grid)) # this should be .0427939
#
# # chi2_test 2:
# a_c2_list = data_to_num_list(a_count_2)
# b_c2_list = data_to_num_list(b_count_2)
# c2_observed_grid = [a_c2_list, b_c2_list]
# print(chi2_value(c2_observed_grid)) # this should be 33.86444
# print(perform_chi2_homogeneity_test(c2_observed_grid)) # this should be 0.0000
# # Again, why do you think this is? Take a peek at a_count_2 and b_count_2 in abtesting_test.py :)
#
# # chi2_test 3:
# a_c3_list = data_to_num_list(a_count_3)
# b_c3_list = data_to_num_list(b_count_3)
# c3_observed_grid = [a_c3_list, b_c3_list]
# print(chi2_value(c3_observed_grid)) # this should be .3119402
# print(perform_chi2_homogeneity_test(c3_observed_grid)) # this should be .57649202
