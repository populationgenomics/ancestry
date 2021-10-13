"""Test Michael's script to understand hail batch"""

import hailtop.batch as hb
import pandas as pd
from google.cloud import storage

OUTPUT_BUCKET = 'cpg-tob-wgs-test'
DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:a3fcb20ebe6cbf75794f736956d1f7b927317cc9-hail-0.2.73.devc6f6f09cec08'  # noqa: E501; pylint: disable=line-too-long

# Define variable OUTSIDE the function,
# hail actually allows you to include this
x = 5  # pylint: disable=invalid-name


def data_frame_with_10_and_var(y):
    """create and return a dataframe"""
    return pd.DataFrame([y, 10, x])


def df_plus_5(df):
    """perform operation and return"""
    return df + 5


def write_dataframe(df, output_path):
    """write df using google.cloud"""
    client = storage.Client()
    bucket = client.get_bucket(OUTPUT_BUCKET)
    bucket.blob(output_path).upload_from_string(df.to_csv(), 'text/csv')


backend = hb.ServiceBackend()
b = hb.Batch(name='run-python job', backend=backend, default_python_image=DRIVER_IMAGE)

# Create the first python_job to create the dataframe
j1 = b.new_python_job(name='create-dataframe')

# this is a pointer to a dataframe, it's actually a Resource (not a writable-resource)
# In the background, hail batch serializes this to a file using `dill`.
df1: hb.resource.PythonResult = j1.call(data_frame_with_10_and_var, 5)

j2 = b.new_python_job('perform-operation-on-dataframe')
# We pass df1 to this function call, or the PythonResult resource,
# Hail batch deserializes this in the backgorund into a dataframe using `dill`.
df2: hb.resource.PythonResult = j2.call(df_plus_5, df1)

# If we want to export it using the df methods to a proper file,
# we should do it in a python_job, eg:
j3 = b.new_python_job('write-dataframe')
j3.call(write_dataframe, df2, 'out-dataframe.csv')
# ,0
# 0,10
# 1,15
# 2,10

# Otherwise hail batch offers us a few alternatives:
#   - `PythonResult.as_str`: `str(object)`
#   - `PythonResult.as_json`: `json.dumps(object)`
#   - `PythonResult.as_repr`: `repr(object)`
b.write_output(df1.as_str(), f'gs://{OUTPUT_BUCKET}/output-df1-str.txt')
#     0
# 0   5
# 1  10
# 2   5

b.write_output(df2.as_repr(), f'gs://{OUTPUT_BUCKET}/output-df2-repr.txt')
#     0
# 0  10
# 1  15
# 2  10


# this fails, as df isn't serializable
# b.write_output(df2.as_json(), "gs://cpg-michael-hail-dev/output-df2.json")

b.run()
