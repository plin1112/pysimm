"""
The pysimm.utils module contains utility classes and functions used throughout the pysimm package.
This includes the Item and ItemContainer classes for managing simulation components,
and the PysimmError exception class for error handling.
"""
# ******************************************************************************
# pysimm.utils module
# ******************************************************************************
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2016 Michael E. Fortunato, Coray M. Colina
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from __future__ import print_function
from collections.abc import Sequence
import itertools
from pprint import pprint
from typing import List, Dict, Optional, Any, Tuple, Union, Iterator

from pysimm import debug_print, warning_print


class PysimmError(Exception):
    """Custom exception class for pySIMM errors."""
    pass


class Container(object):
    """
    An arbitrary container object that returns None if trying to access an attribute that does not exist.
    """
    def __getattr__(self, name: str) -> None:
        return None


class Item(object):
    """
    A generic class for objects in a simulation, such as particles, bonds, and angles.
    It allows for arbitrary attributes to be set via keyword arguments.
    """
    def __init__(self, **kwargs: Any) -> None:
        """
        Initializes an Item object, setting attributes from keyword arguments.

        Args:
            **kwargs: Arbitrary keyword arguments to be set as attributes.
        """
        for k, v in kwargs.items():
            if v == 'TRUE':
                v = True
            elif v == 'FALSE':
                v = False
            setattr(self, k, v)

    def __getattr__(self, name: str) -> Any:
        """
        Returns None if an attribute is not found.

        Args:
            name: The name of the attribute.

        Returns:
            The attribute value or None if not found.
        """
        return None

    def copy(self) -> 'Item':
        """
        Creates a shallow copy of the item.

        Returns:
            A new Item object with the same attributes.
        """
        i = type(self)()
        for k, v in vars(self).items():
            setattr(i, k, v)
        if hasattr(i, 'tag') and i.tag:
            del i.tag
        return i

    def set(self, **kwargs: Any) -> None:
        """
        Sets attributes from keyword arguments.

        Args:
            **kwargs: Arbitrary keyword arguments to be set as attributes.
        """
        for k, v in kwargs.items():
            setattr(self, k, v)


class ItemContainer(Sequence):
    """
    A container for organizing Item objects, such as particles, bonds, and angles.
    It uses a dictionary for storage, with integer tags as keys.
    """
    def __init__(self, _dict: Optional[Dict[int, Item]] = None, **kwargs: Any) -> None:
        """
        Initializes an ItemContainer.

        Args:
            _dict: An optional dictionary of items to initialize the container with.
            **kwargs: Arbitrary keyword arguments to be set as attributes.
        """
        self._dict: Dict[int, Item] = _dict or {}
        self.count: int = len(self._dict) if _dict else 0
        self._name_map: Dict[str, List[Item]] = {}
        if _dict:
            for item in self._dict.values():
                self._add_to_name_map(item)
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __len__(self) -> int:
        return len(self._dict)

    def __iter__(self) -> Iterator[Item]:
        for v in list(self._dict.values()):
            yield v

    def __getitem__(self, slice_: Union[int, slice]) -> Union[Item, List[Item], None]:
        if not isinstance(slice_, slice) and isinstance(slice_, int):
            if slice_ >= 0:
                return self._dict.get(slice_)
            else:
                return self._dict.get(self.count + slice_ + 1)
        elif isinstance(slice_, slice):
            data: List[Item] = []
            start, stop, step = slice_.indices(len(self._dict))
            for i in range(start+1, stop+1, step):
                item = self._dict.get(i)
                if item:
                    data.append(item)
            return data
        else:
            return None

    def _add_to_name_map(self, item: Item) -> None:
        if hasattr(item, 'name') and item.name:
            if item.name not in self._name_map:
                self._name_map[item.name] = []
            self._name_map[item.name].append(item)
        if hasattr(item, 'rname') and item.rname and item.rname != item.name:
            if item.rname not in self._name_map:
                self._name_map[item.rname] = []
            self._name_map[item.rname].append(item)

    def _remove_from_name_map(self, item: Item) -> None:
        if hasattr(item, 'name') and item.name and item.name in self._name_map:
            try:
                self._name_map[item.name].remove(item)
                if not self._name_map[item.name]:
                    del self._name_map[item.name]
            except ValueError:
                pass  # item was not in list
        if hasattr(item, 'rname') and item.rname and item.rname != item.name and item.rname in self._name_map:
            try:
                self._name_map[item.rname].remove(item)
                if not self._name_map[item.rname]:
                    del self._name_map[item.rname]
            except ValueError:
                pass  # item was not in list

    def add(self, _item: Item) -> Optional[Item]:
        """
        Adds an item to the container.

        Args:
            _item: The item to add.

        Returns:
            The added item, or None if the addition failed.
        """
        if hasattr(_item, 'tag') and _item.tag is not None and self._dict.get(_item.tag) is None:
            self._dict[_item.tag] = _item
            self.count += 1
        elif (not hasattr(_item, 'tag') or _item.tag is None) and not self._dict.get(self.count+1):
            self.count += 1
            self._dict[self.count] = _item
            _item.tag = self.count
        else:
            debug_print('cannot add at index %s' % (self.count+1))
            return None
        self._add_to_name_map(_item)
        return _item

    def get(self, *args: Any, **kwargs: Any) -> List[Item]:
        """
        Retrieves items from the container based on various criteria.

        Args:
            *args: A single positional argument can be used to specify the name of the item to retrieve.
            **kwargs: Keyword arguments for more specific queries (e.g., tags, startswith).

        Returns:
            A list of items that match the query.
        """
        name: Optional[str] = None
        tags: Optional[List[int]] = kwargs.get('tags')
        startswith: Optional[str] = kwargs.get('startswith')
        first: bool = kwargs.get('first', False)
        order: bool = kwargs.get('order', False)
        item_wildcard: str = kwargs.get('item_wildcard', 'X')
        query_wildcard: Optional[str] = kwargs.get('query_wildcard')
        improper_type: bool = kwargs.get('improper_type', False)
        if len(args) == 1:
            name = args[0]
            if name == 'all':
                return list(self._dict.values())

        if name and not query_wildcard and (not item_wildcard or item_wildcard not in name):
            results = self._name_map.get(name, [])
            if order:
                found = [item for item in results if hasattr(item, 'name') and item.name == name]
            else:
                found = results

            if first and found:
                return [found[0]]
            return sorted(found, key=lambda x: x.name.count('X') if hasattr(x, 'name') and x.name else 0)

        if (len(args) != 1 and kwargs.get('startswith') is None or
                kwargs.get('startswith') and kwargs.get('tags')):
            warning_print('improper usage of system.ItemContainer.get')
            return []
        found: List[Item] = []
        for item in self:
            item_name = getattr(item, 'name', None)
            item_rname = getattr(item, 'rname', None)
            if name is not None and (item_name == name or
                                     (not order and item_rname == name)):
                found.append(item)
                if first:
                    return found

            elif name:
                if item_name and compare(name, item_name, query_wildcard=query_wildcard,
                           item_wildcard=item_wildcard,
                           improper_type=improper_type, order=order):
                    found.append(item)
                    if first:
                        return found
            elif startswith:
                if item_name and item_name.startswith(startswith):
                    found.append(item)
                elif (not item_name and hasattr(item, 'type') and hasattr(item.type, 'name') and
                        item.type.name.startswith(startswith)):
                    found.append(item)
                if first:
                    return found
            elif tags and isinstance(tags, list):
                if hasattr(item, 'tag') and item.tag in tags:
                    found.append(item)
                    if first:
                        return found
        return sorted(found, key=lambda x: x.name.count('X') if hasattr(x, 'name') and x.name else 0)

    def remove(self, index: Union[int, str], update: bool = True) -> None:
        """
        Removes an item from the container.

        Args:
            index: The tag of the item to remove, or 'all' to remove all items.
            update: If True, re-indexes the container to maintain contiguous tags.
        """
        if index == 'all':
            self._dict.clear()
            self._name_map.clear()
            self.count = 0
            return

        if isinstance(index, str):
            raise PysimmError('string index for remove is not supported')

        item_to_remove = self._dict.get(index)
        if item_to_remove:
            self._remove_from_name_map(item_to_remove)
            if not update or index == self.count:
                del self._dict[index]
            else:
                for i in range(index, self.count):
                    item_to_move = self._dict.pop(i + 1, None)
                    if item_to_move:
                        self._dict[i] = item_to_move
                        self._dict[i].tag = i
                    else:
                        raise PysimmError(f'Cannot find item with tag {i+1} to move during re-indexing.')

                if self.count in self._dict:
                    del self._dict[self.count]
            self.count -= 1


def compare(query: str, item: str, query_wildcard: Optional[str] = None, item_wildcard: str = 'X', order: bool = False,
            improper_type: bool = False) -> bool:
    """
    Compares two strings with optional wildcard support.

    Args:
        query: The query string.
        item: The item string to compare against.
        query_wildcard: An optional wildcard character for the query string.
        item_wildcard: The wildcard character for the item string.
        order: If True, the order of comma-separated parts matters.
        improper_type: If True, handles improper type comparisons.

    Returns:
        True if the strings match, False otherwise.
    """
    query_parts = query.split(',')
    item_parts = item.split(',')

    if len(query_parts) != len(item_parts):
        debug_print('cannot compare names with different lengths: %s %s'
                    % (query, item))
        return False

    match = []
    if improper_type:
        if query_parts[0] != item_parts[0]:
            return False
        for q_perm in itertools.permutations(query_parts[1:]):
            for i_perm in itertools.permutations(item_parts[1:]):
                if compare(','.join(q_perm), ','.join(i_perm)):
                    return True
        return False

    for q, i in zip(query_parts, item_parts):
        if q == i:
            match.append(True)
        elif (item_wildcard and i == item_wildcard) or \
                (query_wildcard and q == query_wildcard):
            match.append(True)
        else:
            match.append(False)
    if all(match):
        return True

    if not order:
        match = []
        rquery = list(reversed(query_parts))
        for q, i in zip(rquery, item_parts):
            if q == i:
                match.append(True)
            elif (item_wildcard and i == item_wildcard) or \
                    (query_wildcard and q == query_wildcard):
                match.append(True)
            else:
                match.append(False)
        if all(match):
            return True

    return False
