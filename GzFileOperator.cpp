#include "GzFileOperator.h"

#include <fstream>
#include <string>

#include <zlib/zlib.h>

GzFileOperator::GzFileOperator(const std::string &gz_file_path, const OpenMode open_mode)
{
	m_gz_file = nullptr;
	m_open_mode = open_mode;
	gzFile tmp_gz_file = gzopen(gz_file_path.c_str(), open_mode == RB ? "rb" : "wb");
	if (tmp_gz_file != nullptr)
	{
		m_gz_file = tmp_gz_file;
	}

}
GzFileOperator::~GzFileOperator()
{
	if (isLoaded())
	{
		gzclose(m_gz_file);
		m_gz_file = nullptr;
	}
}
bool GzFileOperator::reopen(const std::string &gz_file_path, const OpenMode open_mode)
{
	if (m_gz_file != nullptr)
	{
		gzclose(m_gz_file);
		m_gz_file = nullptr;
	}
	m_open_mode = open_mode;
	gzFile tmp_gz_file = gzopen(gz_file_path.c_str(), open_mode == RB ? "rb" : "wb");
	if (tmp_gz_file != nullptr)
	{
		m_gz_file = tmp_gz_file;
	}
}
bool GzFileOperator::isLoaded()
{
	return m_gz_file != nullptr;
}
bool GzFileOperator::uncompressTo(const std::string &target_file_path)
{
	// 局部化成员属性
	gzFile gz_file = m_gz_file;
	OpenMode open_mode = m_open_mode;
	if (!isLoaded())
	{
		// 若本文件(源)未加载，则返回
		return false;
	}
	if (open_mode == WB)
	{
		// 若本文件(源)打开模式为写二进制，则返回(应当为读二进制)
		return false;
	}
	// 加载目标文件
	std::ofstream target_file{ target_file_path, std::ios_base::out | std::ios_base::binary };
	if (!target_file.is_open())
	{
		// 若目标文件打开失败，则返回
		return false;
	}
	char *buffer = new char[kMAX_BUFFER_LEN];
	char buffer_debug[1000];
	unsigned int len = kMAX_BUFFER_LEN;
	int actual_read_len;
	do 
	{
		actual_read_len = gzread(gz_file, static_cast<voidp>(buffer_debug), len); // 从压缩文件中读取未压缩数据
		if (actual_read_len != -1)
		{
			target_file.write(buffer_debug, actual_read_len); // 将未压缩数据写入目标文件中
		}
		else
		{
			// 数据错误，返回
			return false;
		}
	} while (actual_read_len==len);

	delete[] buffer;
	target_file.close();
	return true;
}
bool GzFileOperator::compressFrom(const std::string &source_file_path)
{
	return false;
	return true;
}