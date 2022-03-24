#pragma once

#include <string>

#include <zlib/zlib.h>

class GzFileOperator
{
public:
	typedef enum OpenMode
	{
		RB = 0,
		WB = 1
	} OpenMode;

public:
	GzFileOperator(const std::string &gz_file_path, const OpenMode open_mode=RB);
	~GzFileOperator();

public:
	bool reopen(const std::string &gz_file_path, const OpenMode open_mode = RB);
	bool isLoaded();
	bool uncompressTo(const std::string &target_file_path);
	// TODO: 
	bool compressFrom(const std::string &source_file_path);
	

private:
	// Ω˚”√øΩ±¥ππ‘Ï
	GzFileOperator(const GzFileOperator &gz_file_operator);

private:
	const unsigned int kMAX_BUFFER_LEN = 1000;

	gzFile m_gz_file;
	OpenMode m_open_mode;
};

