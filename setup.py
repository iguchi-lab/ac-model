from setuptools import setup, find_packages

setup(
    name="ac_model",  # パッケージ名
    version="0.1.0",  # バージョン
    packages=find_packages(include=["ac_model", "ac_model.*"]),
    install_requires=[
        "requests",
        "numpy",
        "pandas",
        "archenv",  # GitHub リンクを直接書かない
    ],
    dependency_links=[
        "https://github.com/iguchi-lab/archenv.git#egg=archenv",
    ],
    include_package_data=True,
    zip_safe=False,
)
