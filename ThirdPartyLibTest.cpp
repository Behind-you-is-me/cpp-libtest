/*
Copyright 2020 Esri

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of
the License at http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

For additional information, contact:
Environmental Systems Research Institute, Inc.
Attn: Contracts Dept
380 New York Street
Redlands, California, USA 92373
email: contracts@esri.com
*/

//#include "i3s/i3s_writer.h"
//#include "utils/utl_png.h"
//#include "utils/utl_geom.h"
//#include "utils/utl_i3s_resource_defines.h"
//#include <iostream>
//#include <vector>
//#include <algorithm>
//#include <cmath>
//#include <cstdlib>
//#include <string>
//
//#include <filesystem>
//namespace stdfs = std::filesystem;
//
//using i3slib::utl::Vec2i;
//using i3slib::utl::Vec2f;
//using i3slib::utl::Vec2d;
//using i3slib::utl::Vec3f;
//using i3slib::utl::Vec3d;
//
//namespace
//{
//
//	constexpr double c_pi = 3.14159265358979323846; // M_PI;
//	constexpr double c_wgs84_equatorial_radius = 6378137.0; //semi-major 
//	constexpr double c_wgs84_polar_radius = 6356752.314245; //semi-minor
//	constexpr double c_degrees_per_meter = 180.0 / (c_wgs84_polar_radius * c_pi);
//
//	struct CS_transformation
//	{
//		DECL_PTR(CS_transformation);
//		virtual bool transform(Vec3d* points, size_t count) = 0;
//	};
//
//	class ENU_to_WGS_transformation : public CS_transformation
//	{
//	public:
//
//		DECL_PTR(ENU_to_WGS_transformation);
//
//		explicit ENU_to_WGS_transformation(const Vec2d& origin) :
//			origin_(origin),
//			lon_degrees_per_meter_(180.0 / (c_wgs84_equatorial_radius * c_pi * std::cos(origin.y * c_pi / 180.0)))
//		{}
//
//		virtual bool transform(Vec3d* points, size_t count) override
//		{
//			std::transform(points, points + count, points, [this](const Vec3d& p)
//			{
//				return Vec3d(origin_.x + p.x * lon_degrees_per_meter_, origin_.y + p.y * c_degrees_per_meter, p.z);
//			});
//
//			return true;
//		}
//
//	private:
//
//		const Vec2d origin_;
//		const double lon_degrees_per_meter_;
//	};
//
//	double screen_size_to_area(double pixels)
//	{
//		constexpr double c_pi_over_4 = c_pi * 0.25;
//		return pixels * pixels * c_pi_over_4;
//	}
//
//	i3slib::i3s::Layer_writer::Var create_writer(const stdfs::path& slpk_path)
//	{
//		i3slib::i3s::Ctx_properties ctx_props;
//		//i3slib::i3s::set_geom_compression(ctx_props.geom_encoding_support, i3slib::i3s::Geometry_compression::Draco, true);
//		//i3slib::i3s::set_gpu_compression(ctx_props.gpu_tex_encoding_support, i3slib::i3s::GPU_texture_compression::ETC_2, true);
//		auto writer_context = i3slib::i3s::create_i3s_writer_context(ctx_props);
//
//		i3slib::i3s::Layer_meta meta;
//		meta.type = i3slib::i3s::Layer_type::Mesh_IM;
//		meta.name = slpk_path.stem().u8string();
//		meta.desc = "Generated with raster2slpk";
//		meta.sr.wkid = 4326;
//		meta.uid = meta.name;
//		meta.normal_reference_frame = i3slib::i3s::Normal_reference_frame::Not_set;
//
//		std::unique_ptr<i3slib::i3s::Layer_writer> writer(
//			i3slib::i3s::create_mesh_layer_builder(writer_context, slpk_path));
//
//		if (writer)
//			writer->set_layer_meta(meta);
//		return writer;
//	}
//
//	template<typename T>
//	size_t remove_alpha_channel(T* data, size_t size)
//	{
//		I3S_ASSERT(size % 4 == 0);
//		size_t out = 3;
//		for (size_t in = 4; in < size; in += 4, out += 3)
//		{
//			data[out] = data[in];
//			data[out + 1] = data[in + 1];
//			data[out + 2] = data[in + 2];
//		}
//
//		I3S_ASSERT(out * 4 == size * 3);
//		return out;
//	}
//
//	bool load_color_data(const stdfs::path& path, int& size, std::vector<char>& data)
//	{
//		data.clear();
//
//		int h;
//		if (!i3slib::utl::read_png_from_file(path, &size, &h, &data))
//		{
//			std::cout << "Failed to load color bitmap from file." << std::endl;
//			return false;
//		}
//
//		if (size != h)
//		{
//			std::cout << "Color bitmap must have equal width and height." << std::endl;
//			return false;
//		}
//
//		if ((size & (size - 1)) != 0)
//		{
//			std::cout << "Color bitmap size must be a power of 2." << std::endl;
//			return false;
//		}
//
//		// The current png reader implementation always produces RGBA output for color images.
//		// We have to strip alpha channel here.
//		I3S_ASSERT(data.size() == size * size * 4);
//		data.resize(remove_alpha_channel(data.data(), data.size()));
//		I3S_ASSERT(data.size() == size * size * 3);
//		return true;
//	}
//
//	bool load_elevation_data(const stdfs::path& path, int size, std::vector<double>& data, double unit)
//	{
//		data.clear();
//
//		int w, h;
//		std::vector<char> grayscale_data;
//		if (!i3slib::utl::read_png_from_file(path, &w, &h, &grayscale_data))
//		{
//			std::cout << "Failed to load elevation image from file." << std::endl;
//			return false;
//		}
//
//		if (w != size + 1 || h != size + 1)
//		{
//			std::cout << "Invalid elevation image dimensions." << std::endl;
//			return false;
//		}
//
//		if (grayscale_data.size() != w * w * 2)
//		{
//			std::cout << "Elevation image is not a 16-bit grayscale." << std::endl;
//			return false;
//		}
//
//		data.reserve(grayscale_data.size() / 2);
//		auto p = reinterpret_cast<const unsigned char*>(grayscale_data.data());
//		auto end = p + grayscale_data.size();
//		for (; p != end; p += 2)
//		{
//			// PNG files store 16-bit pixels in network byte order (that is big-endian, most significant bits
//			// first). See http://www.libpng.org/pub/png/libpng-manual.txt
//			// However, implementation of read_png_from_file() sets png_set_invert_mono() and png_set_swap()
//			// modes for libpng, and we have to deal with this here.
//
//			const auto e = 0xffffu - (p[0] + (p[1] << 8));
//			data.push_back(e * unit);
//		}
//
//		return true;
//	}
//
//	void build_downsampled_textures(int size, const char* data, int min_size, std::vector<std::vector<char>>& textures)
//	{
//		while (size > min_size)
//		{
//			I3S_ASSERT((size % 2) == 0);
//			const auto s = size / 2;
//			const auto stride = size * 3;
//
//			std::vector<char> texture;
//			texture.reserve(s * s * 3);
//
//			const unsigned char *p = reinterpret_cast<const unsigned char*>(data);
//			for (int y = 0; y < s; y++)
//			{
//				auto p1 = p + stride;
//				for (int x = 0; x < s; x++, p += 6, p1 += 6)
//				{
//					texture.push_back(static_cast<char>(std::round((p[0] + p[3] + p1[0] + p1[3]) * 0.25)));
//					texture.push_back(static_cast<char>(std::round((p[1] + p[4] + p1[1] + p1[4]) * 0.25)));
//					texture.push_back(static_cast<char>(std::round((p[2] + p[5] + p1[2] + p1[5]) * 0.25)));
//				}
//				p = p1;
//			}
//
//			textures.emplace_back(std::move(texture));
//			data = textures.back().data();
//			size = s;
//		}
//	}
//
//	void build_downsampled_grids(int size, const double* data, int min_size, std::vector<std::vector<double>>& grids)
//	{
//		while (size > min_size)
//		{
//			I3S_ASSERT((size % 2) == 0);
//			const auto s = size / 2;
//			I3S_ASSERT((s % 2) == 0);
//
//			std::vector<double> grid;
//			grid.reserve((s + 1) * (s + 1));
//
//			auto p = data;
//			auto pn = p + size + 1;
//
//			// Corner vertex of the first row.
//			grid.push_back((p[0] + p[1] + pn[0]) / 3.0);
//			p += 2;
//			pn += 2;
//
//			// Internal vertices of the first row.
//			for (int x = 1; x < s; x++, p += 2, pn += 2)
//				grid.push_back((p[-1] + p[0] + p[1] + pn[0]) * 0.25);
//
//			// Corner vertex of the first row.
//			grid.push_back((p[-1] + p[0] + pn[0]) / 3.0);
//
//			//
//			auto pp = p + 1;
//			p = pn + 1;
//			pn += (size + 2);
//
//			for (int y = 1; y < s; y++)
//			{
//				// First vertex of the row.
//				grid.push_back((pp[0] + p[0] + p[1] + pn[0]) * 0.25);
//				pp += 2;
//				p += 2;
//				pn += 2;
//
//				// Internal vertices of the row.
//				for (int x = 1; x < s; x++, pp += 2, p += 2, pn += 2)
//					grid.push_back((pp[0] + p[-1] + p[0] + p[1] + pn[0]) * 0.2);
//
//				// Last vertex of the row.
//				grid.push_back((pp[0] + p[-1] + p[0] + pn[0]) * 0.25);
//
//				pp = p + 1;
//				p = pn + 1;
//				pn += (size + 2);
//			}
//
//			// Corner vertex of the last row.
//			grid.push_back((pp[0] + p[0] + p[1]) / 3.0);
//			pp += 2;
//			p += 2;
//
//			// Internal vertices of the last row.
//			for (int x = 1; x < s; x++, pp += 2, p += 2, pn += 2)
//				grid.push_back((pp[0] + p[-1] + p[0] + p[1]) * 0.25);
//
//			// Corner vertex of the last row.
//			grid.push_back((pp[0] + p[-1] + p[0]) / 3.0);
//
//			grids.emplace_back(std::move(grid));
//			data = grids.back().data();
//			size = s;
//		}
//	}
//
//	//double get_elevation(const double* grid_data, int grid_size, int x, int y)
//	//{
//	//  return grid_data[y * (grid_size + 1) + x];
//	//};
//
//	void add_quad(
//		std::vector<Vec3d>& verts,
//		std::vector<Vec2f>& uvs,
//		const Vec3d& v0,
//		const Vec2f& uv0,
//		const Vec3d& v1,
//		const Vec2f& uv1,
//		const Vec3d& v2,
//		const Vec2f& uv2,
//		const Vec3d& v3,
//		const Vec2f& uv3)
//	{
//		verts.push_back(v0);
//		verts.push_back(v1);
//		verts.push_back(v2);
//		verts.push_back(v2);
//		verts.push_back(v3);
//		verts.push_back(v0);
//
//		uvs.push_back(uv0);
//		uvs.push_back(uv1);
//		uvs.push_back(uv2);
//		uvs.push_back(uv2);
//		uvs.push_back(uv3);
//		uvs.push_back(uv0);
//	}
//
//	// Extracts the fragment of the elevation grid starting at _start_ and spanning _size_ cells 
//	// (i.e, size + 1 nodes) and extracts the fragment of color_data starting at texture_start
//	// with texture_size * texture_size pixels.
//
//	bool build_mesh(
//		const i3slib::i3s::Layer_writer& writer,
//		const Vec2d& cell_size,
//		int cell_factor,
//		const double* grid_data,
//		int grid_size,
//		const char* color_data,
//		int color_size,
//		const Vec2i& start,
//		int size,
//		const Vec2i& texture_start,
//		int texture_size,
//		i3slib::i3s::Mesh_data& mesh,
//		CS_transformation* transformation = nullptr)
//	{
//		I3S_ASSERT(start.x + size <= grid_size);
//		I3S_ASSERT(start.y + size <= grid_size);
//		I3S_ASSERT(texture_start.x + texture_size <= color_size);
//		I3S_ASSERT(texture_start.y + texture_size <= color_size);
//
//		const auto end = start + Vec2i(size, size);
//		const double h = cell_factor * cell_size.x;
//
//		const auto get_elevation = [&grid_data, grid_size](int x, int y)
//		{
//			return grid_data[y * (grid_size + 1) + x];
//		};
//
//		std::vector<Vec3d> verts;
//		std::vector<Vec2f> uvs;
//		verts.reserve(static_cast<size_t>(size + 3) * (size + 3) * 6);
//		uvs.reserve(verts.capacity());
//
//		if (start.y != 0)
//		{
//			// Add a row of skirt quads.
//			Vec2f uv(0, 0);
//			Vec2f uv_next(0, 0);
//			for (int dx = 0; dx < size; dx++, uv.x = uv_next.x)
//			{
//				auto ind_x = start.x + dx;
//				const Vec3d vtx0((ind_x * cell_factor) * cell_size.x, -(start.y * cell_factor) * cell_size.y, get_elevation(ind_x, start.y));
//				const Vec3d vtx1(((ind_x + 1) * cell_factor) * cell_size.x, vtx0.y, get_elevation(ind_x + 1, start.y));
//				const Vec3d vtx2(vtx1.x, vtx1.y, vtx1.z - h);
//				const Vec3d vtx3(vtx0.x, vtx0.y, vtx0.z - h);
//				uv_next.x = static_cast<float>(dx + 1) / size;
//				add_quad(verts, uvs, vtx0, uv, vtx1, uv_next, vtx2, uv_next, vtx3, uv);
//			}
//		}
//
//		for (int dy = 0; dy < size; dy++)
//		{
//			const auto ind_y = start.y + dy;
//			const auto ind_y1 = ind_y + 1;
//			const auto y = -(ind_y * cell_factor) * cell_size.y;
//			const auto y1 = -(ind_y1 * cell_factor) * cell_size.y;
//			const auto v = static_cast<float>(dy) / size;
//			const auto v1 = static_cast<float>(dy + 1) / size;
//
//			if (start.x != 0)
//			{
//				// Add skirt quad at the beginning of the row.
//				const Vec3d vtx0((start.x * cell_factor) * cell_size.x, y1, get_elevation(start.x, ind_y1));
//				const Vec3d vtx1(vtx0.x, y, get_elevation(start.x, ind_y));
//				const Vec3d vtx2(vtx1.x, vtx1.y, vtx1.z - h);
//				const Vec3d vtx3(vtx0.x, vtx0.y, vtx0.z - h);
//				add_quad(verts, uvs, vtx0, { 0, v1 }, vtx1, { 0, v }, vtx2, { 0, v }, vtx3, { 0, v1 });
//			}
//
//			for (int dx = 0; dx < size; dx++)
//			{
//				const auto ind_x = start.x + dx;
//				const auto ind_x1 = ind_x + 1;
//				const Vec3d v00((ind_x * cell_factor) * cell_size.x, y, get_elevation(ind_x, ind_y));
//				const Vec3d v11((ind_x1 * cell_factor) * cell_size.x, y1, get_elevation(ind_x1, ind_y1));
//				const Vec3d v10(v11.x, v00.y, get_elevation(ind_x1, ind_y));
//				const Vec3d v01(v00.x, v11.y, get_elevation(ind_x, ind_y1));
//				const Vec2f uv00(static_cast<float>(dx) / size, v);
//				const Vec2f uv11(static_cast<float>(dx + 1) / size, v1);
//				const Vec2f uv01(uv00.x, v1);
//				const Vec2f uv10(uv11.x, v);
//				add_quad(verts, uvs, v10, uv10, v00, uv00, v01, uv01, v11, uv11);
//			}
//
//			if (end.x != grid_size)
//			{
//				// Add skirt quad at the end of the row.
//				const Vec3d vtx0((end.x * cell_factor) * cell_size.x, y, get_elevation(end.x, ind_y));
//				const Vec3d vtx1(vtx0.x, y1, get_elevation(end.x, ind_y1));
//				const Vec3d vtx2(vtx0.x, y1, vtx1.z - h);
//				const Vec3d vtx3(vtx0.x, y, vtx0.z - h);
//
//				const Vec2f uv0(1.0f, v);
//				const Vec2f uv1(1.0f, v1);
//
//				// The 0.9999 hack is due to a probable bug in Pro SLPK renderer.
//				// If a triangle has all its u texture coordinates equal to 1.0, it is textured incorrectly.
//				// Web Scene Viewer does not have this issue.
//				add_quad(verts, uvs, vtx0, uv0, vtx1, uv1, vtx2, /*uv1*/{ 0.9999f, v1 }, vtx3, /*uv0*/{ 0.9999f, v });
//			}
//		}
//
//		if (end.y != grid_size)
//		{
//			// Add a row of skirt quads.
//			Vec2f uv(0, 1.0f);
//			Vec2f uv_next(0, 1.0f);
//			for (int dx = 0; dx < size; dx++, uv.x = uv_next.x)
//			{
//				auto ind_x = start.x + dx;
//				const Vec3d vtx0(((ind_x + 1) * cell_factor) * cell_size.x, -(end.y * cell_factor) * cell_size.y, get_elevation(ind_x + 1, end.y));
//				const Vec3d vtx1((ind_x * cell_factor) * cell_size.x, vtx0.y, get_elevation(ind_x, end.y));
//				const Vec3d vtx2(vtx1.x, vtx1.y, vtx1.z - h);
//				const Vec3d vtx3(vtx0.x, vtx0.y, vtx0.z - h);
//				uv_next.x = static_cast<float>(dx + 1) / size;
//
//				// The 0.9999 hack is due to a probable bug in Pro SLPK renderer.
//				// If a triangle has all its v texture coordinates equal to 1.0, it is textured incorrectly.
//				// Web Scene Viewer does not have this issue.
//				add_quad(verts, uvs, vtx0, uv_next, vtx1, uv, vtx2, /*uv*/{ uv.x, 0.9999f }, vtx3, /*uv_next*/{ uv_next.x, 0.9999f });
//			}
//		}
//
//		if (transformation)
//			transformation->transform(verts.data(), verts.size());
//
//		std::vector<char> texture(static_cast<size_t>(texture_size) * texture_size * 3);
//		const auto bytes_per_row = texture_size * 3;
//		const auto stride = color_size * 3;
//		auto row = color_data + (texture_start.y * color_size + texture_start.x) * 3;
//		auto rgb_out = texture.data();
//		for (size_t y = 0; y < texture_size; y++, row += stride)
//			rgb_out = std::copy(row, row + bytes_per_row, rgb_out);
//
//		i3slib::i3s::Simple_raw_mesh raw_mesh;
//		raw_mesh.vertex_count = static_cast<int>(verts.size());
//		raw_mesh.abs_xyz = verts.data();
//		raw_mesh.uv = uvs.data();
//		if (!i3slib::i3s::create_texture_from_image(texture_size, texture_size, 3, texture.data(), raw_mesh.img))
//			return false;
//
//		return writer.create_mesh_from_raw(raw_mesh, mesh) == IDS_I3S_OK;
//	}
//
//	bool process(
//		i3slib::i3s::Layer_writer& writer,
//		const int input_size,
//		const Vec2d& cell_size,
//		const std::vector<std::vector<double>>& grids,
//		const std::vector<std::vector<char>>& textures,
//		int node_tris_size,
//		int node_texture_size,
//		int depth,
//		int grid_size,
//		const Vec2i& start,
//		int texture_size,
//		const Vec2i& texture_start,
//		i3slib::i3s::Node_id& node_id,
//		CS_transformation* transformation = nullptr)
//	{
//		std::vector<i3slib::i3s::Node_id> node_ids;
//
//		if (depth + 1 < textures.size())
//		{
//			for (int i : { 0, 1 })
//			{
//				for (int j : { 0, 1 })
//				{
//					const auto status = process(
//						writer, input_size, cell_size, grids, textures, node_tris_size, node_texture_size, depth + 1,
//						grid_size * 2, 2 * start + node_tris_size * Vec2i(j, i),
//						texture_size * 2, 2 * texture_start + node_texture_size * Vec2i(j, i),
//						node_id, transformation);
//
//					if (!status)
//						return false;
//
//					node_ids.push_back(node_id++);
//				}
//			}
//		}
//		else if (depth + 1 < grids.size())
//		{
//			for (int i : { 0, 1 })
//			{
//				for (int j : { 0, 1 })
//				{
//					const auto status = process(
//						writer, input_size, cell_size, grids, textures, node_tris_size, node_texture_size / 2, depth + 1,
//						grid_size * 2, 2 * start + node_tris_size * Vec2i(j, i),
//						texture_size, texture_start + node_texture_size / 2 * Vec2i(j, i),
//						node_id, transformation);
//
//					if (!status)
//						return false;
//
//					node_ids.push_back(node_id++);
//				}
//			}
//		}
//
//		//
//		i3slib::i3s::Simple_node_data node_data;
//		node_data.node_depth = depth + 1;
//		node_data.children = std::move(node_ids);
//		node_data.lod_threshold = screen_size_to_area(500);
//
//		const std::vector<char>& texture = depth < textures.size() ? textures[depth] : textures.back();
//
//		const auto status = build_mesh(
//			writer, cell_size, input_size / grid_size,
//			grids[depth].data(), grid_size, texture.data(), texture_size,
//			start, node_tris_size, texture_start, node_texture_size, node_data.mesh, transformation);
//
//		if (!status)
//			return false;
//
//		return writer.create_node(node_data, node_id) == IDS_I3S_OK;
//	}
//
//}
//
//int main(/*int argc, char* argv[]*/)
//{
//	int argc = 7;
//	const char* argv[7] = { "path", R"(D:\PersonalData\DataSet\ps_height_1k.png)", R"(D:\PersonalData\DataSet\ps_texture_1k.png)", R"(.\普吉特海湾_1k.slpk)", "160", "160", "0.1" };
//	if (argc != 7)
//	{
//		std::cout << "Usage:" << std::endl
//			<< "raster2slpk <elevation_png> <color_png> <output_slpk_file> <x_step> <y_step> <z_unit>" << std::endl;
//
//		return 1;
//	}
//
//	const stdfs::path elevation_file_path(argv[1]);
//	const stdfs::path color_file_path(argv[2]);
//	const stdfs::path slpk_file_path(argv[3]);
//
//	const Vec2d cell_size(std::stod(argv[4]), std::stod(argv[5]));
//	const double elevation_unit = std::stod(argv[6]);
//
//	//
//	int size;
//	std::vector<char> color_data;
//	if (!load_color_data(color_file_path, size, color_data))
//		return 1;
//
//	std::vector<double> elevation_data;
//	if (!load_elevation_data(elevation_file_path, size, elevation_data, elevation_unit))
//		return 1;
//
//	std::vector<std::vector<char>> textures;
//	textures.emplace_back(std::move(color_data));
//	build_downsampled_textures(size, textures.front().data(), 128, textures);
//	std::reverse(std::begin(textures), std::end(textures));
//
//	std::vector<std::vector<double>> grids;
//	grids.emplace_back(std::move(elevation_data));
//	build_downsampled_grids(size, grids.front().data(), 32, grids);
//	std::reverse(std::begin(grids), std::end(grids));
//
//	//
//	auto writer = create_writer(slpk_file_path);
//	if (!writer)
//		return 1;
//
//	ENU_to_WGS_transformation transformation({ -123.4583943, 47.6204856 });
//	i3slib::i3s::Node_id node_id = 0;
//
//	if (!process(*writer, size, cell_size, grids, textures, 32, 128, 0, 32, { 0, 0 }, 128, { 0, 0 }, node_id, &transformation))
//		return 1;
//
//	// Add a root node on top of everything.
//	i3slib::i3s::Simple_node_data node_data;
//	node_data.node_depth = 0;
//	node_data.children.push_back(node_id++);
//	if (writer->create_node(node_data, node_id) != IDS_I3S_OK)
//		return 1;
//
//	if (writer->save() != IDS_I3S_OK)
//		return 1;
//
//	return 0;
//}
//























// ThirdPartyLibTest.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

#include <json/json.h>
#include <zlib/zlib.h>
#include <zip.h>
#include <sqlapi/SQLAPI.h>
#include <i3s-lib/i3s/i3s_enums.h>

#include "GzFileOperator.h"
#include "AccessDBCast.h"

void JsonTest();
void ZLibTest();
void ZLibRead();
void ZLibWrite();
void ZipTest();
bool extractSLPKTo(const std::string &slpk_file_path, const std::string &target_dir);
bool extractGZTo(const std::string &slpk_file_path, const std::string &target_dir);
std::string exeSQL(const std::string &sql,const SAConnection *sa_con);
void SQLTest();
void TestGzFileOperator();
void DBUpdateTest();

int main(int argn, char** args)
{
	/*JsonTest();
	ZLibTest();
	ZLibRead();
	ZLibWrite();*/
	//ZipTest();
	//DBUpdateTest();
	std::cout << "argn: " << argn << std::endl;
	for (int arg_i = 0; arg_i < argn; ++arg_i)
	{
		std::cout << arg_i << ": " << *(args + arg_i) << std::endl;
	}

	SAConnection sa_conn;
	try
	{
		sa_conn.Connect("LH_Access_SQL_Test", "", "", SA_ODBC_Client);
		std::string help_info = "\n1. 输入 exit + Enter 退出程序.\n"
			"2. 整行输入\n\n";
		std::cout << help_info;
		std::string sql;
		while (true)
		{
			std::cout << ">> ";
			std::getline(std::cin, sql);
			if (sql == "exit")
			{
				std::cout << sql << std::endl;
				break;
			}
			std::cout << "result: " << exeSQL(sql, &sa_conn) << "\n";

		}
	}
	catch (const SAException&e)
	{
		std::string err_txt = e.ErrText().GetMultiByteChars();
		std::cout << err_txt << std::endl;
		std::cout << "exit.\n";
	}

	//std::string path = R"(E:\Desktop\tmp\Buildings_NewYork\nodes\1\attributes\f_2\0.bin)";
	//std::ifstream fin;
	//fin.open(path, std::ios_base::in | std::ios_base::binary);
	//UINT32 i, count, oid=0;
	//double SOURCE_ID;
	//int bytes = sizeof(SOURCE_ID);
	//fin.read((char*)(&count), 4);
	//bytes = sizeof(oid);
	////fin.read((char*)(&SOURCE_ID), 4);
	//std::map< double, UINT32> cc;
	//for (i = 0; i < count; ++i)
	//{
	//	if (fin.eof())
	//	{
	//		break;
	//	}
	//	fin.read((char*)(&oid), bytes);
	//	++cc[oid];
	//}
	//
	//bool end_of_file = fin.eof();
	//for (; !fin.eof();)
	//{
	//	fin.read((char*)(&oid), bytes);
	//	++i;
	//}
	//fin.close();
	//std::cout << "end";
	//SQLTest();
	//TestGzFileOperator();
	return 0;
}
void JsonTest()
{
	Json::Value root;
	root["name"] = "Youle";
	root["age"] = 21;
	root["salary"] = 429.67;
	root["hobby"] = "Food, Game, and another things that is interested.";
	Json::Value task;
	task.append(4);
	task.append(7.990);
	task.append("real life.");
	root["task"] = task;

	std::string json_string = Json::StyledWriter().write(root);
	std::ofstream fout("tmp/test.json");
	fout << json_string;
	fout.close();
}
std::string getResultTxtOfZLib(const int status_code)
{
	switch (status_code)
	{
	case Z_OK:
		return "Z_OK";
	case Z_MEM_ERROR:
		return "Z_MEM_ERROR";
	case Z_BUF_ERROR:
		return "Z_BUF_ERROR";
	case Z_DATA_ERROR:
		return "Z_DATA_ERROR";
	default:
		return "Unknown_Error";
	}
}
void ZLibTest()
{
	const uLong MAX_LEN = 1000;

	std::string gz_file_path = "tmp/0_0_1.bin.dds.gz";
	Bytef dest[MAX_LEN]= "1234";
	uLong dest_len = MAX_LEN;
	std::string data = "jfdsapghnasdp giupdsahjg pidsipafghipuas gipsiujhgpiuaqjp";
	int status_code = compress(dest, &dest_len, (const Bytef*)data.c_str(), data.size());
	std::string result_txt = getResultTxtOfZLib(status_code);
	
	Bytef uncompress_data[MAX_LEN] = "1234";
	uLong uncompress_data_len = 10;
	status_code = uncompress(uncompress_data, &uncompress_data_len, dest, dest_len);
	result_txt = getResultTxtOfZLib(status_code);
}
void ZLibRead()
{
	gzFile file = gzopen("tmp/0_0_1.bin.dds.gz", "rb");
	if (file == nullptr)
	{
		return;
	}
	const uLong MAX_LEN = 100;
	char data[MAX_LEN];
	voidp buf = data;
	unsigned len = MAX_LEN;
	int readed_size = gzread(file, buf, len);
	std::ofstream fout{ "tmp/pic.dds", std::ios_base::out | std::ios_base::binary };
	while (readed_size == len)
	{
		fout.write(data, readed_size);
		readed_size = gzread(file, buf, len);
	}
	gzclose(file);
}
void ZLibWrite()
{
	// 二进制方式读文本。gzputs, ok；gzwrite, err
	std::ifstream fin{ "tmp/test.json", std::ios_base::in | std::ios_base::binary };
	const int MAX_LEN = 100;
	char data[MAX_LEN];
	gzFile gz_file = gzopen("tmp/test_copied.json.gz", "wb");
	unsigned int data_len = 0;
	int writed_len = 0;
	gzFile gz_file_err = gzopen("tmp/test_copied_err.json.gz", "wb");
	while (!fin.eof())
	{
		fin.read(data, MAX_LEN-1);
		data_len = fin.gcount();
		data[data_len] = '\0';
		writed_len = gzwrite(gz_file_err, (voidpc)data,data_len);
		writed_len = gzputs(gz_file, data);
	}
	int status_code = gzclose(gz_file);
	std::string status_str = getResultTxtOfZLib(status_code);
	status_code = gzclose(gz_file_err);
	status_str = getResultTxtOfZLib(status_code);
	fin.close();
	// 二进制方式读二进制。写二进制ok；不可用写文本的方式(因为空字符的数据文本结束标识)
	fin.open("tmp/0_0_1.bin.dds", std::ios_base::in | std::ios_base::binary);
	data[MAX_LEN];
	gz_file = gzopen("tmp/0_0_1_copied.dds.gz", "wb");
	data_len = 0;
	writed_len = 0;
	while (!fin.eof())
	{
		fin.read(data, MAX_LEN - 1);
		data_len = fin.gcount();
		writed_len = gzwrite(gz_file, (voidpc)data,data_len);
		//writed_len = gzputs(gz_file, data);
	}
	status_code = gzclose(gz_file);
	status_str = getResultTxtOfZLib(status_code);
}
void ZipTest()
{
	// 解压缩文件
	// 1. 打开档案
	std::string zip_path = "tmp/testF5.slpk";
	int error_code;
	zip_t *archive_ptr = zip_open(zip_path.c_str(), ZIP_RDONLY, &error_code);
	if (archive_ptr == nullptr)
	{
		return;
	}
	// 2. 获取档案中的文件个数(名字)
	int file_count = zip_get_num_files(archive_ptr);
	std::string file_name = "";
	zip_file_t *zip_file_ptr = nullptr;
	const int MAX_BUFFER_SIZE = 100;
	char Buffer[MAX_BUFFER_SIZE];
	void *buffer_ptr = Buffer;
	zip_int64_t read_byte_count = 0;
	std::ofstream fout;
	for (int file_index = 0; file_index < file_count; ++file_index)
	{
		file_name = zip_get_name(archive_ptr, file_index, ZIP_FL_ENC_GUESS);
		// 保险一点，确认文件在档案中存在
		zip_int64_t real_index = zip_name_locate(archive_ptr, file_name.c_str(), ZIP_FL_ENC_GUESS);
		if (real_index == -1)
		{
			continue;
		}
		// 3. 根据文件名或文件索引读取文件
		zip_file_ptr = zip_fopen(archive_ptr, file_name.c_str(), ZIP_FL_COMPRESSED);
		//zip_file_prt = zip_fopen_index(archive_ptr, file_index, ZIP_FL_COMPRESSED);
		std::filesystem::path dir_path{ "tmp/extract/" + file_name };
		fout.open("tmp/extract/" + file_name, std::ios_base::out | std::ios_base::binary);
		if (*dir_path.string().rbegin() == '/')
		{
			std::filesystem::create_directories(dir_path);
			continue;
		}
		do
		{
			read_byte_count = zip_fread(zip_file_ptr, buffer_ptr, MAX_BUFFER_SIZE);
			// 4. 将数据写入文件中
			fout.write(Buffer, read_byte_count);
		} while (read_byte_count == MAX_BUFFER_SIZE);
		fout.close();
	}

	// 压缩文件
}
bool extractSLPKTo(const std::string &slpk_file_path, const std::string &target_dir)
{
	std::filesystem::path target_path{ target_dir };
	if (!std::filesystem::exists(target_path))
	{
		std::filesystem::create_directories(target_dir);
	}
	// 解压缩文件
	// 1. 打开档案
	int error_code;
	zip_t *archive_ptr = zip_open(slpk_file_path.c_str(), ZIP_RDONLY, &error_code);
	if (archive_ptr == nullptr)
	{
		return false;
	}
	// 2. 获取档案中的文件个数(名字)
	int file_count = zip_get_num_files(archive_ptr);
	std::string file_name = "";
	zip_file_t *zip_file_ptr = nullptr;
	const int MAX_BUFFER_SIZE = 8*1024;
	char Buffer[MAX_BUFFER_SIZE];
	void *buffer_ptr = Buffer;
	zip_int64_t read_byte_count = 0;
	std::ofstream fout;
	for (int file_index = 0; file_index < file_count; ++file_index)
	{
		file_name = zip_get_name(archive_ptr, file_index, ZIP_FL_ENC_GUESS);
		// 保险一点，确认文件在档案中存在
		zip_int64_t real_index = zip_name_locate(archive_ptr, file_name.c_str(), ZIP_FL_ENC_GUESS);
		if (real_index == -1)
		{
			continue;
		}
		// 3. 根据文件名或文件索引读取文件
		zip_file_ptr = zip_fopen(archive_ptr, file_name.c_str(), ZIP_FL_COMPRESSED);
		//zip_file_prt = zip_fopen_index(archive_ptr, file_index, ZIP_FL_COMPRESSED);
		std::filesystem::path dir_path{ target_dir +"/" + file_name };
		
		if (*dir_path.string().rbegin() == '/')
		{
			std::filesystem::create_directories(dir_path);
			continue;
		}
		else if(!std::filesystem::exists(dir_path.parent_path()))
		{
			std::filesystem::create_directories(dir_path.parent_path());
		}
		fout.open(dir_path.string(), std::ios_base::out | std::ios_base::binary);
		do
		{
			read_byte_count = zip_fread(zip_file_ptr, buffer_ptr, MAX_BUFFER_SIZE);
			// 4. 将数据写入文件中
			fout.write(Buffer, read_byte_count);
		} while (read_byte_count == MAX_BUFFER_SIZE);
		fout.close();
		// 将gz文件解压缩一份用以查看
		if (dir_path.string().substr(dir_path.string().size()-3) == ".gz")
		{
			bool success = extractGZTo(dir_path.string(), dir_path.string().substr(0, dir_path.string().size() - dir_path.filename().string().size() - 1));
		}
	}
	return true;
}
bool extractGZTo(const std::string &slpk_file_path, const std::string &target_dir)
{
	gzFile file = gzopen(slpk_file_path.c_str(), "rb");
	if (file == nullptr)
	{
		return false;
	}
	const uLong MAX_LEN = 8 * 1024;
	char data[MAX_LEN];
	voidp buf = data;
	unsigned len = MAX_LEN;
	int readed_size = 0;
	std::filesystem::path file_path{ slpk_file_path };
	std::string file_name = file_path.filename().string();
	std::string target_file_name_path = target_dir + "/" + file_name.substr(0, file_name.size() - 3);
	std::ofstream fout{ target_file_name_path, std::ios_base::out | std::ios_base::binary };
	do
	{
		readed_size = gzread(file, buf, len);
		fout.write(data, readed_size);
	} while (readed_size == len);
	gzclose(file);
	return true;
}
std::string exeSQL(const std::string &sql, const SAConnection *sa_con_ptr)
{
	SACommand sa_cmd;
	try
	{
		sa_cmd.setConnection(const_cast<SAConnection *>(sa_con_ptr));
		sa_cmd.setCommandText(sql.c_str());
		sa_cmd.Execute();
		return "execute ok.";
	}
	catch (const SAException& e)
	{
		std::string err_txt = e.ErrText().GetMultiByteChars();
		return err_txt;
	}
	
}
void SQLTest()
{
	SAConnection sa_con;
	try
	{
		sa_con.Connect("LH_Access_sql_test", "", "", SA_ODBC_Client);
	}
	catch (const SAException&e)
	{
		std::string err_txt = e.ErrText().GetMultiByteChars();
		std::cout << err_txt << std::endl;
		return;
	}
	std::string sql = "insert into text_longt(t_Lt) values(:1);";
	std::string value = R"(灰黄色，岩石风化剧烈，组织结构基本破坏，仅局部可辨。矿物组成以石英颗粒和长石、云母 为主，长石、云母等易风化矿物已基本风化成次生粘土矿物。岩芯呈坚硬土状，该岩具浸水软化崩解，力学强度降低的工程特性。此层与上部的残积土层无明显界限，是按实测标贯锤击数30≤N＜50击确定的，该层岩体结构类型为散体状结构，岩石为极软岩、岩体极破碎， 岩石质量指标RQD为0，其等级属"极差的"，岩体基本质量等级为Ⅴ级)";
	SACommand sa_cmd;
	try
	{
		sa_cmd.setConnection(const_cast<SAConnection *>(&sa_con));
		sa_cmd.setCommandText(sql.c_str());
		sa_cmd.Param("1").setAsString() = SAString(value.c_str());
		sa_cmd.Execute();
	}
	catch (const SAException& e)
	{
		std::string err_txt = e.ErrText().GetMultiByteChars();
	}
	
	/*while (true)
	{
		std::getline(std::cin, sql);
		if (sql == "exit")
		{
			break;
		}
		std::cout << exeSQL(sql, &sa_con) << std::endl;
	}*/

}
void TestGzFileOperator()
{
	GzFileOperator gz_operator{ "tmp/testF5.slpk" };
	gz_operator.uncompressTo("tmp/testF5_test");

}
void DBUpdateTest()
{
	SAConnection rf_db_conn;
	try
	{
		rf_db_conn.Connect("LH_Access_target", "", "", SA_ODBC_Client);
	}
	catch (const SAException& e)
	{
		std::string err_txt = e.ErrText().GetMultiByteChars();
	}
	AccessDBCast access_db_cast(R"(E:\Desktop\tmp\xm-zd-db\target)", &rf_db_conn);
	access_db_cast.castToTarget(R"(E:\Desktop\tmp\xm-zd-db\厦门数据库空.mdb)");
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
