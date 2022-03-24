#include "AccessDBCast.h"

#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <sqlapi/SQLAPI.h>

void initReferenceInfo(std::map< std::string, std::vector<std::pair< std::string, SADataType_t>> >& ref_info, const SAConnection* rf_db_ptr)
{
	ref_info.clear();
	if (rf_db_ptr == nullptr)
	{
		return;
	}
	SACommand sa_cmd;
	std::string sql;
	try
	{
		sa_cmd.setConnection(const_cast<SAConnection*>(rf_db_ptr));
		sql = "select name from MSysObjects where (type=1 or type=5)  and flags =0";
		sa_cmd.setCommandText(sql.c_str());
		sa_cmd.Execute();
		
		std::vector<std::string> all_table_name;
		while (sa_cmd.FetchNext())
		{
			all_table_name.push_back(sa_cmd.Field(1).asString().GetMultiByteChars());
		}
		
		for (std::string &table_name : all_table_name)
		{
			sql = "select * from " + table_name + " where false";
			sa_cmd.setCommandText(sql.c_str());
			sa_cmd.Execute();
			std::vector<std::pair<std::string, SADataType_t>> &field_type = ref_info[table_name];
			int field_count = sa_cmd.FieldCount();
			for (int i = 0; i < field_count; ++i)
			{
				field_type.push_back({ sa_cmd.Field(i + 1).Name().GetMultiByteChars(), sa_cmd.Field(i + 1).FieldType() });
			}
		}

	}
	catch (const SAException&e)
	{
		std::string err_txt = e.ErrText().GetMultiByteChars();
	}
}
bool copyFileTo(const std::string &source_file_path, const std::string &target_file_path, const unsigned int byte_chunk = 2*1024*1024/* 2MB */)
{
	std::ifstream source_stream;
	std::ofstream target_stream;
	source_stream.open(source_file_path, std::ios_base::in | std::ios_base::binary);
	target_stream.open(target_file_path, std::ios_base::out | std::ios_base::binary);
 
	std::streamsize byte_read = 0;
	char *buffer_ptr = new char[byte_chunk];
	while (!source_stream.eof())
	{
		source_stream.read(buffer_ptr, byte_chunk);
		byte_read = source_stream.gcount();
		target_stream.write(buffer_ptr, byte_read);
	}
	delete[] buffer_ptr;
	target_stream.close();
	source_stream.close();
	return true;
}

AccessDBCast::AccessDBCast(const std::string &target_dir, const SAConnection* rf_db_ptr)
{
	m_target_dir = target_dir;
	initReferenceInfo(m_reference_info, rf_db_ptr);

}
bool AccessDBCast::castToTarget(const std::string &source_access_file, std::string *err_txt) // TODO: DB_CONN_INFO
{
	auto setErrTxt = [&err_txt](const std::string &msg) {
		err_txt != nullptr ? (*err_txt = msg, 1) : 0;
	};
	std::string target_file_path(m_target_dir);
	// �õ�Ŀ���ļ�·��
	std::string file_name = "target";
	char separator = '/';
	size_t pos = source_access_file.find_last_of(separator);
	if (pos == std::string::npos)
	{
		separator = '\\';
		pos = source_access_file.find_last_of(separator);
	}
	if (pos == std::string::npos)
	{
		setErrTxt("�ļ�·�����Ϸ����޷������ļ���");
		return false; // �����ļ�·��
	}
	file_name = source_access_file.substr(pos + 1);
	if ((*target_file_path.rbegin()) != separator)
	{
		target_file_path.push_back(separator);
	}
	target_file_path += file_name;
	// ����Դ�ļ���Ŀ���ļ�
	copyFileTo(source_access_file, target_file_path);
	// �޸�Ŀ���ļ�
	auto isFieldCompatible = [](std::vector<std::pair< std::string, SADataType_t>>::const_iterator &found_it,
		const std::pair<std::string, SADataType_t> &field_info,
		const std::vector<std::pair< std::string, SADataType_t>> &ref_info)->bool {
		// �ڲο���Ϣ�в��Ҵ��ֶεĲο���Ϣ
		for (found_it = ref_info.begin(); found_it != ref_info.end() && found_it->first != field_info.first; ++found_it);
		
		if (found_it == ref_info.end() || found_it->second == field_info.second)
		{
			return true; // û�вο���Ϣ���ֶβ��������������ݴ����ֶ�������ͬ��һ�����ݡ�
		}
		/*
			// TBD: �����ݻ��ǲ�һ���� -- LH.2022.03.09 15:07 wait to decided.
				����˵��(int, double, string)��
					int --> double     ok.
					int --> string     ok.
					double --> int     no.
					string --> int     no.
					string --> double  no.
		*/
		return false;
	};
	auto toString = [](const SADataType_t data_type)->std::string {
		switch (data_type)
		{
		case SADataType_t::SA_dtBool: return "Yes/No";
		case SADataType_t::SA_dtShort:
		case SADataType_t::SA_dtUShort: return "Integer";
		case SADataType_t::SA_dtLong:
		case SADataType_t::SA_dtULong: return "Long";
		case SADataType_t::SA_dtInt64:
		case SADataType_t::SA_dtUInt64:

		case SADataType_t::SA_dtDouble:
		case SADataType_t::SA_dtNumeric: return "Double";
		case SADataType_t::SA_dtDateTime: return "Date/Time";
		case SADataType_t::SA_dtString: return "Text";
		case SADataType_t::SA_dtLongChar: return "Memo";
		default:
			break;
		}
		return "";
	};
	SAConnection sa_conn;
	try
	{
		sa_conn.Connect("LH_Access_target", "", "", SA_ODBC_Client);
		SACommand select_cmd(&sa_conn);
		//SACommand update_cmd(&sa_conn);
		std::string sql = "";
		
		std::vector< std::pair<std::string, SADataType_t> > field_info_list;
		for (auto table_it = m_reference_info.begin(); table_it != m_reference_info.end(); ++table_it)
		{ // �Բο���Ϣ�е�ÿһ�ű�
			sql = "select * from " + table_it->first + " where false";
			select_cmd.setCommandText(sql.c_str());
			try
			{
				select_cmd.Execute();
			}
			catch (const SAException&e)
			{ // ��ѯ��ʧ�ܣ������Ǳ����ڣ��Ǿ�������������ѯ���ű�
				std::string db_err_txt = e.ErrText().GetMultiByteChars();
				continue;
			}
			
			int field_count = select_cmd.FieldCount();
			std::vector<std::pair< std::string, SADataType_t>>::const_iterator ref_found_it;
			field_info_list.clear();
			for (int field_i = 1; field_i <= field_count; ++field_i)
			{ // �ȴ�˱���ֶ���Ϣ
				field_info_list.push_back({ select_cmd.Field(field_i).Name().GetMultiByteChars(), select_cmd.Field(field_i).FieldType() });
			}
			for (auto field_info_it = field_info_list.begin(); field_info_it != field_info_list.end(); ++field_info_it)
			{ // �Դ˱��ÿһ���ֶ�
				// �ж�����ֶ�����ͬ�Ĳο���Ϣ�е��ֶε���������ֶ������Ƿ����
				if (!isFieldCompatible(ref_found_it, *field_info_it, table_it->second))
				{ // �����ݣ�����ǿ���޸Ĵ��ֶ�����Ϊ�ο���Ϣ�е�����
					sql = "ALTER TABLE " + table_it->first + " ALTER COLUMN " + field_info_it->first + " " + toString(ref_found_it->second);
					select_cmd.setCommandText(sql.c_str());
					select_cmd.Execute();
				}
			}
		}
	}
	catch (const SAException& e)
	{
		std::string db_err_txt = e.ErrText().GetMultiByteChars();
		setErrTxt(db_err_txt);
		return false;
	}
	
	return true;
}
