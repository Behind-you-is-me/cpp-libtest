#pragma once

#include <string>
#include <vector>
#include <map>

#include <sqlapi/SQLAPI.h>

class AccessDBCast
{
public:
	AccessDBCast() = delete;
	AccessDBCast(const std::string &target_dir, const SAConnection* rf_db_ptr);

	bool castToTarget(const std::string &source_access_file, std::string *err_txt = nullptr); // TODO: DB_CONN_INFO

private:
	std::string m_target_dir;
	std::map< std::string, std::vector<std::pair< std::string, SADataType_t>> > m_reference_info;
};

