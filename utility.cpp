#include "utility.h"

std::ostream& operator<<(std::ostream& out, const Index_update value) {
	switch (value)
	{
	case Index_update::withlimit: return out << "withlimit";
	case Index_update::withlimit_base_opt: return out << "withlimit_base_opt";
	case Index_update::withlimit_dfs: return out << "withlimit_dfs";
	case Index_update::withlimit_dfs_parallel: return out << "withlimit_dfs_parallel";
	case Index_update::withlimit_parallel: return out << "withlimit_parallel";
	case Index_update::withoutlimit: return out << "withoutlimit";
	};
	return out << static_cast<std::uint16_t>(value);
}