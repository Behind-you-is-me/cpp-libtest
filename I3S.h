#pragma once

#include <string>
#include <filesystem>

#include <zip.h>

class gmeModel;

class I3S
{
public:
	I3S(const std::string &file_path);
	~I3S();

public:
	bool saveFromGmeModel(gmeModel *gme_model);
	gmeModel * toGmeModel();

	bool extractTo(const std::string &folder_path);

private:

private:


};

