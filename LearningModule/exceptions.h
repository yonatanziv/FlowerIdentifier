#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <exception>
#include <string>
using namespace std;

class ContourFailBinaryBackgroundException: public exception
{
	virtual const char* what() const throw()
	{
		return "CONTOUR_FAIL_BINARY_BACKGROUND if there are more white pixels than black \
			   in the area masked by background mask (fail of checkBinaryBackground)\n";
	}
};

class ContourFailCenterPointException: public exception
{
	virtual const char* what() const throw()
	{
		return "Failed to find contour\n";
	}
};

class NoOuterMaxPointsFoundException: public exception
{
	virtual const char* what() const throw()
	{
		return "NO_OUTER_MAX_POINTS_FOUND_ERROR\n";
	}
};

class NoOuterMinPointsFoundException: public exception
{
	virtual const char* what() const throw()
	{
		return "NO_OUTER_MIN_POINTS_FOUND_ERROR\n";
	}
};

class MinMaxPointsNumberException: public exception
{
	virtual const char* what() const throw()
	{
		return "ERROR: min or max points is too little\n";
	}
};

class MedianAngleException: public exception
{
	virtual const char* what() const throw()
	{
		return "ERROR: Bad median angle.\n";
	}
};

class MaxInnerContourException: public exception
{
	virtual const char* what() const throw()
	{
		return "ERROR: Bad median angle.\n";
	}
};

class NoInnerMinMaxPointsException: public exception
{
	virtual const char* what() const throw()
	{
		return "ERROR: Bad median angle.\n";
	}
};

class PathNotFoundException: public exception
{
	virtual const char* what() const throw()
	{
		return "ERROR: path does not exists!\n";

	}
};

class EmptyImageException: public exception
{
	virtual const char* what() const throw()
	{
		return "ERROR: Flower image couldnt be opened\n";

	}
};

#endif