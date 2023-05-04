#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
################################################################################
This script defines the data structure for stack and queue.

Created on 29 September 2021.
################################################################################
"""
__author__ = 'ZLiang'

class StackUnderflow(ValueError):
    # stack under flow (visit empty stack)
    pass

class Stack:
    def __init__(self):
        # Use list to store elements
        self._elements = []

    def is_empty(self):
        return self._elements == []

    def push(self, element):
        self._elements.append(element)

    def pop(self):
        if self._elements == []:
            raise StackUnderflow("In the Stack.pop()")
        return self._elements.pop()

    def top(self):
        if self._elements == []:
            raise StackUnderflow("In the Stack.top()")
        return self._elements[-1]
