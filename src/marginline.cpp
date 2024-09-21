#include "marginline.h"
#include <cassert>
#include <cmath>
#include "curvature_info.h"


namespace
{
	double PerpendicularDistance(
		const Eigen::Vector3d& point,
		const Eigen::Vector3d& line_start,
		const Eigen::Vector3d& line_end
	) {
		Eigen::Vector3d line_vec = line_end - line_start;
		Eigen::Vector3d point_vec = point - line_start;
		auto line_length = line_vec.norm();
		
		if (line_length < 1e-4)	// TODO: tolerance should be defined as a constant
		{
			return point_vec.norm();
		}

		auto t = point_vec.dot(line_vec) / (line_length * line_length);
		Eigen::Vector3d projection = line_start + t * line_vec;
		return (projection - point).norm();
	}


	void RDP(
		const VectorArray& V,
		const std::vector<int>& polyline_vertices,
		double epsilon,
		int start,
		int end,
		std::vector<int>& result
	) {
		if (start >= end)
		{
			return;
		}

		// find the point with the maximum distance from the line connecting the start and end points
		auto max_distance = 0.0;
		auto index = start;
		for (auto i = start + 1; i < end; ++i)
		{
			auto distance = PerpendicularDistance(
				V.row(polyline_vertices[i]),
				V.row(polyline_vertices[start]),
				V.row(polyline_vertices[end]));
			if (distance > max_distance)
			{
				max_distance = distance;
				index = i;
			}
		}

		// keep the point if the maximum distance exceeds the threshold
		if (max_distance > epsilon)
		{
			RDP(V, polyline_vertices, epsilon, start, index, result);
			result.push_back(polyline_vertices[index]);
			RDP(V, polyline_vertices, epsilon, index, end, result);
		}
	}
}


void CreateMarginline(
	const VectorArray& V,
	const IndicesArray& F,
	const std::vector<std::vector<int>>& adjacency_list,
	const CurvatureInfo& curvature_info,
	std::vector<int>& marginline,
	std::set<int>& visited)
{
	static const size_t MAX_NUM_TRAVERSAL = 10000;
	static const int64_t NUM_HOPS = 10;

	if (marginline.empty())
	{
		return;
	}

	visited.clear();
	visited.insert(marginline.begin(), marginline.end());

	auto start = marginline.back();
	for (size_t i = 0; i < MAX_NUM_TRAVERSAL; ++i)
	{
		if (marginline.size() > 1)
		{
			if (marginline.front() == marginline.back())
			{
				break;
			}
		}

		auto seed = marginline.back();
		auto neighbors = adjacency_list[seed];

		const Eigen::Vector3d max_curvature_direction = curvature_info.principal_directions1.row(seed);
		const Eigen::Vector3d min_curvature_direction = curvature_info.principal_directions2.row(seed);

		{
			using IndexAndMeanCurvature = std::tuple<int, double>;
			std::vector< IndexAndMeanCurvature> candidates;
			for (size_t j = 0; j < neighbors.size(); ++j)
			{
				auto& neigbor = neighbors[j];
				auto found = visited.find(neigbor);
				if (found != visited.end())
				{
					continue;
				}

				Eigen::Vector3d direction = (V.row(neigbor) - V.row(seed)).normalized();
				assert(marginline.size() > 0);
				auto start = std::max(static_cast<int64_t>(0), static_cast<int64_t>(marginline.size()) - NUM_HOPS - 1);
				auto end = static_cast<int64_t>(marginline.size() - 1);
				auto is_opposite_direction = false;
				for (int64_t k = start; k < end; ++k)
				{
					Eigen::Vector3d existing_direction = (V.row(marginline[k + 1]) - V.row(marginline[k])).normalized();
					if (direction.dot(existing_direction) < 0.0)
					{
						is_opposite_direction = true;
						break;
					}
				}
				if (is_opposite_direction)
				{
					continue;
				}

				candidates.push_back(std::make_tuple(neigbor, curvature_info.mean[neigbor]));
			}
			if (!candidates.empty())
			{
				auto max_element = std::max_element(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs)
					{
						return std::get<1>(lhs) < std::get<1>(rhs);
					});
				if (std::get<1>(*max_element) > curvature_info.mean[seed])
				{
					seed = std::get<0>(*max_element);
					marginline.push_back(seed);
					visited.insert(neighbors.begin(), neighbors.end());
					continue;
				}
			}
		}

		{
			using IndexAndAbsCos = std::tuple<int, double>;
			std::vector<IndexAndAbsCos> candidates;
			for (size_t j = 0; j < neighbors.size(); ++j)
			{
				auto& neigbor = neighbors[j];
				auto found = visited.find(neigbor);
				if (found != visited.end())
				{
					continue;
				}

				if (curvature_info.mean[seed] > 0 && curvature_info.mean[neigbor] < 0)
				{
					continue;
				}

				Eigen::Vector3d direction = (V.row(neigbor) - V.row(seed)).normalized();
#if 0
				assert(result.size() > 0);
				auto start = std::max(static_cast<int64_t>(0), static_cast<int64_t>(result.size()) - NUM_HOPS - 1);
				auto end = static_cast<int64_t>(result.size() - 1);
				auto is_opposite_direction = false;
				for (int64_t k = start; k < end; ++k)
				{
					Eigen::Vector3d existing_direction = (V.row(result[k + 1]) - V.row(result[k])).normalized();
					if (direction.dot(existing_direction) < 0.0)
					{
						is_opposite_direction = true;
						break;
					}
				}
				if (is_opposite_direction)
				{
					continue;
				}
#endif

				candidates.push_back(std::make_tuple(neigbor, std::abs(direction.dot(min_curvature_direction))));
			}
			if (candidates.empty())
			{
				break;
			}

			auto next = std::max_element(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs)
				{
					return std::get<1>(lhs) < std::get<1>(rhs);
				});
			seed = std::get<0>(*next);
			marginline.push_back(seed);
			visited.insert(neighbors.begin(), neighbors.end());
		}
	}
}


std::vector<int> DownSampleMarginline(
	const VectorArray& V,
	const std::vector<int>& marginline,
	size_t num_samples,
	double threshold_to_remove_last_point)
{
	if (marginline.size() < 2)
	{
		return marginline;
	}

    std::vector<int> result;
	static const double epsilon = 0.5;	// TODO: epsilon should be defined as a constant
    result.push_back(marginline.front());
    RDP(V, marginline, epsilon, 0, static_cast<int>(marginline.size() - 1), result);
    result.push_back(marginline.back());

    return result;
}