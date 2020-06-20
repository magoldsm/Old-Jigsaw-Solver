#pragma once

struct BivertexArc
{
	Curve	m_Contour;
	Curve	m_Signature;
};

using BATable = Eigen::Matrix<BivertexArc, Eigen::Dynamic, 1>;

class CPiece
{
public:
	CPiece();
	~CPiece();

	void operator=(std::vector<cv::Point2d>& contour);

	Curve			m_Contour;
	Curve			m_Signature;
	double			m_Weight;
	int				m_nActive;
	BATable			m_Arcs;
	Eigen::VectorXi	m_Pt2Arc;

	size_t		size() { return m_Contour.size(); }
};

bool ReadPuzzle(std::vector<CPiece>& Pieces,  const char* filename);

