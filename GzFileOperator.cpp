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
	// �ֲ�����Ա����
	gzFile gz_file = m_gz_file;
	OpenMode open_mode = m_open_mode;
	if (!isLoaded())
	{
		// �����ļ�(Դ)δ���أ��򷵻�
		return false;
	}
	if (open_mode == WB)
	{
		// �����ļ�(Դ)��ģʽΪд�����ƣ��򷵻�(Ӧ��Ϊ��������)
		return false;
	}
	// ����Ŀ���ļ�
	std::ofstream target_file{ target_file_path, std::ios_base::out | std::ios_base::binary };
	if (!target_file.is_open())
	{
		// ��Ŀ���ļ���ʧ�ܣ��򷵻�
		return false;
	}
	char *buffer = new char[kMAX_BUFFER_LEN];
	char buffer_debug[1000];
	unsigned int len = kMAX_BUFFER_LEN;
	int actual_read_len;
	do 
	{
		actual_read_len = gzread(gz_file, static_cast<voidp>(buffer_debug), len); // ��ѹ���ļ��ж�ȡδѹ������
		if (actual_read_len != -1)
		{
			target_file.write(buffer_debug, actual_read_len); // ��δѹ������д��Ŀ���ļ���
		}
		else
		{
			// ���ݴ��󣬷���
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