#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
#####################################################################################################
This script tests the data structure for stack and queue.

Created on 25 February 2022 for the unit test using pytest.
#####################################################################################################
"""
__author__ = 'ZLiang'

import pytest
from spectral_library.stack_queue import Stack,StackUnderflow

def test_is_empty():
    stack = Stack()
    assert stack.is_empty() is True
    stack.push(12365)
    assert stack.is_empty() is not True

def test_push():
    stack = Stack()
    stack.push(12365)
    stack.push(7)
    assert stack.is_empty() is not True

def test_pop():
    stack = Stack()
    stack.push(12365)
    stack.push(7)
    stack.push("(")
    assert stack.pop() == "("
    assert stack.pop() == 7
    assert stack.pop() == 12365
    try:
        stack.pop()
    except StackUnderflow:
        print("Stack Under flow")

def test_top():
    stack = Stack()
    stack.push(12365)
    assert stack.top() == 12365
    stack.push(7)
    assert stack.top() == 7
    stack.push("(")
    assert stack.top() == "("

