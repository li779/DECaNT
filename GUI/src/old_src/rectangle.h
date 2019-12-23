#ifndef rectangle_h
#define rectangle_h

#include <list>
#include <memory>

using namespace std;

struct position
{
	double x;
	double z;
};

class rectangle
{

	double H; //height of rectangle
	double W; //width of rectangle
	int tubeNum; //the tube that the rectangle is associated with
	position pos; //position of bottom left corner
	double area; //area of the rectangle
	


public:
	inline rectangle();
	inline rectangle(double width, double height);
	inline rectangle(double width, double height, int new_tubeNum);
	inline ~rectangle();
	inline void setWidth(double width);
	inline void setHeight(double height);
	inline void setPosition(position new_pos);
	inline void setPosition(double x, double z);
	inline double getWidth();
	inline double getHeight();
	inline position getPosition();
	inline double getArea();
	inline int getTubeNum();
	inline void setTubeNum(int new_tubeNum);
	
};

struct level
{
	list<shared_ptr<rectangle>> rectsOnLevel;
	double lvlHeight;
	double usedWidth;
	int levelNum;
};


inline rectangle::rectangle()
{
}

inline rectangle::rectangle(double width, double height)
{
	W = width;
	H = height;
	pos.x = 0.0;
	pos.z = 0.0;
	area = W*H;

}

inline rectangle::rectangle(double width, double height, int new_tubeNum)
{
	W = width;
	H = height;
	pos.x = 0.0;
	pos.z = 0.0;
	area = W*H;
	tubeNum = new_tubeNum;
}

inline rectangle::~rectangle()
{
}

inline double rectangle::getHeight()
{
	return H;
}

inline position rectangle::getPosition()
{
	return pos;
}

inline double rectangle::getWidth()
{
	return W;
}

inline void rectangle::setHeight(double height)
{
	H = height;
	area = W*H;
}


inline void rectangle::setPosition(double x, double z)
{
	pos.x = x;
	pos.z = z;
}

inline void rectangle::setPosition(position new_pos)
{
	pos = new_pos;
}

inline void rectangle::setWidth(double width)
{
	W = width;
	area = W*H;
}

inline double rectangle::getArea()
{
	return area;
}

inline int rectangle::getTubeNum()
{
	return tubeNum;
}

inline void rectangle::setTubeNum(int new_tubeNum)
{
	tubeNum = new_tubeNum;
}

#endif //rectangle_h