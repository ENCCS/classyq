/*
 * ClassyQ, classical polarizable solvent models.
 * Copyright (C) 2022 Roberto Di Remigio Eikås and contributors.
 *
 * This file is part of ClassyQ.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * For information on the complete list of contributors to the
 * ClassyQ library, see: <https://classyq.readthedocs.io/>
 */

#pragma once

#include <cstdlib>
#include <string>

#include <fmt/format.h>

namespace classyq {
/** Information on allocation size in human-readable units.
 *
 * @tparam T type of allocated data.
 * @param[in] count number of elements allocated.
 * @return string representation of amount of allocated data.
 */
template<typename T>
auto
memory_with_units(size_t count) -> std::string
{
  std::string units[5] = { "B", "KiB", "MiB", "GiB", "TiB" };

  int i = 0;
  auto bytes = sizeof(T) * count;
  auto mem_bytes = static_cast<double>(bytes);

  for (auto i = 0; i < 5 && bytes >= 1024; ++i, bytes /= 1024) {
    mem_bytes = bytes / 1024.0;
  }

  return fmt::format("{:.2f} {:s}", mem_bytes, units[i]);
}
} // namespace classyq
