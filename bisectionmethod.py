# Python program for implementation 
# of Bisection Method for 
# solving equations 


# An example function whose 
# solution is determined using 
# Bisection Method. 
# The function is x^3 - x^2 + 2


#

mannings_n= 0.011 # check
b_width= 1.5 # meters
depth_assumed= 2.0 # meters
side_slope = 0.0
slope=0.02
def func(depth_assumed):

    A_f= (b_width +side_slope * depth_assumed)*depth_assumed  # in m2
    P_peri= b_width + 2*(side_slope**2+1)*depth_assumed  # in meters
    R_radius= float(A_f)/P_peri  #
    Q_cal= (1/mannings_n)*A_f*(R_radius**(0.666))*slope**(0.5)  # in m3/s
    return 16.0-Q_cal

# Prints root of func(x) 
# with error of EPSILON 
def bisection(a,b): 

	if (func(a) * func(b) >= 0): 
		print("You have not assumed right a and b\n") 
		return

	c = a 
	while (abs(b-a) >= 0.01): 

		# Find middle point 
		c = (a+b)/2

		# Check if middle point is root 
		if (func(c) == 0.0): 
			break

		# Decide the side to repeat the steps 
		if (func(c)*func(a) < 0): 
			b = c 
		else: 
			a = c 
			
	print("The value of root is : ","%.4f"%c) 
	
# Driver code 
# Initial values assumed 
a =20.0
b = 0.0
c=15
d=1
store_1=bisection(a, b) 
store_2=bisection(c,d)

# This code is contributed 
# by Anant Agarwal. 
