from excitingtools.dict_utils import container_converter

def test_convert_container():
    data1 = {'a': '5.0', 'b': ['1.0', '10.0'], 'c': {'a1': '1.0'}, 'd': "blu", 'e': [['1', '2.3'], '3'], 'f': ['1', 'a']}
    data2 = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2.3], 3], 'f': [1, 'a']}

    assert container_converter(data1) == data2

def test_convert_container_no_strings():
    data1 = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2.3], 3], 'f': [1, 11.3]}
    data2 = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2.3], 3], 'f': [1, 11.3]}
    
    assert container_converter(data1) == data2

def test_convert_container_data_type():
    data1 = {'a': '5.0', 'b': ['1.0', '10.0'], 'c': {'a1': '1.0'}, 'd': "blu", 'e': [['1', '2.3'], '3']}
    data2 = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2.3], 3]}

    data1 = container_converter(data1)
    
    for elem1, elem2 in zip(data1.values(), data2.values()):
        assert type(elem1) == type(elem2)
    for elem, elem2 in zip(data1['c'].values(), data2['c'].values()): 
        assert type(elem) == type(elem2)
    for elem, elem2 in zip(data1['b'], data2['b']): 
        assert type(elem) == type(elem2)
    for elem, elem2 in zip(data1['e'], data2['e']): 
        assert type(elem) == type(elem2)
    for elem, elem2 in zip(data1['e'][0], data2['e'][0]):
        assert type(elem) == type(elem2)

    
