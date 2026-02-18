import streamlit as st
import pandas as pd
import numpy as np
import altair as alt

# Set random seed for reproducibility
np.random.seed(42)

st.set_page_config(layout="wide")
st.title("Compact Dashboard App")

# Define dataframes first
df1 = pd.DataFrame({
    "Product": [f"Product {i}" for i in range(1, 1001)],
    "Sales": np.random.randint(100, 1000, 1000),
    "Quarter": np.random.choice(["Q1", "Q2", "Q3", "Q4"], 1000),
    "Region": np.random.choice(["North", "South", "East", "West"], 1000),
    "Profit": np.random.randint(20, 200, 1000),
    "Progress": np.random.uniform(0, 1, 1000),
    "Status": np.random.choice(["✅ Complete", "⏳ In Progress", "❌ Failed"], 1000),
    "Growth": np.random.uniform(-0.3, 0.5, 1000).round(2),
    "Product Link": [f"https://example.com/product/{i}" for i in range(1, 1001)],
    "Chart": np.random.randn(1000, 20).tolist()  # Small sparkline data
})

df2 = pd.DataFrame({
    "Customer": [f"Customer {i}" for i in range(1, 1001)],
    "Avatar": [f"https://i.pravatar.cc/150?img={i%70}" for i in range(1, 1001)],
    "Age": np.random.randint(18, 80, 1000),
    "City": np.random.choice(["NY", "LA", "CHI", "HOU", "SF", "MIA", "SEA", "BOS", "DEN", "ATL"], 1000),
    "Loyalty Level": np.random.choice(["🥉 Bronze", "🥈 Silver", "🥇 Gold", "💎 Diamond"], 1000),
    "Total Spend": np.random.randint(100, 2000, 1000),
    "Satisfaction": np.random.uniform(1.0, 5.0, 1000).round(1),
    "Notes": [f"• Last purchase: Electronics\n• Preferred: Online\n• Support tickets: {n}" 
              for n in np.random.randint(0, 5, 1000)]
})

df3 = pd.DataFrame({
    "Item": [f"Item {i}" for i in range(1, 1001)],
    "Image": [f"https://picsum.photos/150/150?random={i}" for i in range(1, 1001)],
    "Stock Status": np.random.choice(["🟢 In Stock", "🟡 Low Stock", "🔴 Out of Stock"], 1000),
    "Stock": np.random.randint(0, 1000, 1000),
    "Price": np.random.uniform(9.99, 99.99, 1000).round(2),
    "Rating": np.random.uniform(1, 5, 1000).round(1),
    "Category": np.random.choice(["Basic", "Standard", "Premium", "Deluxe"], 1000),
    "Tags": [[f"tag{i}" for i in np.random.choice(range(1,6), np.random.randint(1,4))] 
             for _ in range(1000)]
})

# Create filters in sidebar
with st.sidebar:
    # Sales filters
    st.header("Sales Filters")
    quarters = st.multiselect("Quarter", options=df1["Quarter"].unique())
    regions = st.multiselect("Region", options=df1["Region"].unique())
    status = st.multiselect("Status", options=df1["Status"].unique())
    min_sales = st.number_input("Min Sales ($)", value=0, step=100)
    st.write("---")
    
    # Customer filters  
    st.header("Customer Filters")
    cities = st.multiselect("City", df2["City"].unique())
    loyalty_levels = st.multiselect("Loyalty Level", options=df2["Loyalty Level"].unique())
    age_range = st.slider("Age Range", 
                         min_value=int(df2["Age"].min()),
                         max_value=int(df2["Age"].max()),
                         value=(int(df2["Age"].min()), int(df2["Age"].max())))
    min_satisfaction = st.slider("Min Satisfaction", 1.0, 5.0, 1.0, 0.5)
    st.write("---")
    
    # Inventory filters
    st.header("Inventory Filters")
    price_range = st.slider("Price Range ($)",
                           min_value=float(df3["Price"].min()),
                           max_value=float(df3["Price"].max()),
                           value=(float(df3["Price"].min()), float(df3["Price"].max())))
    categories = st.multiselect("Category", options=df3["Category"].unique())
    stock_status = st.multiselect("Stock Status", options=df3["Stock Status"].unique())

# Apply filters
if quarters:
    df1 = df1[df1["Quarter"].isin(quarters)]
if regions:
    df1 = df1[df1["Region"].isin(regions)]
if status:
    df1 = df1[df1["Status"].isin(status)]
df1 = df1[df1["Sales"] >= min_sales]
    
if cities:
    df2 = df2[df2["City"].isin(cities)]
if loyalty_levels:
    df2 = df2[df2["Loyalty Level"].isin(loyalty_levels)]
df2 = df2[(df2["Age"] >= age_range[0]) & (df2["Age"] <= age_range[1])]
df2 = df2[df2["Satisfaction"] >= min_satisfaction]

if categories:
    df3 = df3[df3["Category"].isin(categories)]
if stock_status:
    df3 = df3[df3["Stock Status"].isin(stock_status)]
df3 = df3[(df3["Price"] >= price_range[0]) & (df3["Price"] <= price_range[1])]

col1, col2, col3 = st.beta_columns(3)

with col1:
    st.write("Sales Data")
    st.dataframe(
        df1[["Product", "Progress", "Chart", "Sales", "Quarter", "Region", "Profit", "Status", "Growth", "Product Link"]],
        # column_config={
        #     "Product": "Product Name",
        #     "Sales": st.column_config.NumberColumn("Sales ($)", format="$%d"),
        #     "Progress": st.column_config.ProgressColumn("Progress", min_value=0, max_value=1),
        #     "Status": "Status",
        #     "Growth": st.column_config.NumberColumn("Growth", format="%+.1f%%"),
        #     "Product Link": st.column_config.LinkColumn("Link"),
        #     "Chart": st.column_config.LineChartColumn("Trend (30d)", width="medium"),
        # },
        # use_container_width=True,
        height=500,
    )

with col2:
    st.write("Customer Data")
    st.dataframe(
        df2,
        # column_config={
        #     "Avatar": st.column_config.ImageColumn("Avatar", width="small"),
        #     "Total Spend": st.column_config.NumberColumn("Total Spend", format="$%d"),
        #     "Last Visit": st.column_config.DateColumn("Last Visit"),
        #     "Satisfaction": st.column_config.NumberColumn("Rating", format="%.1f⭐"),
        #     "Notes": st.column_config.TextColumn("Notes", width="medium"),
        #     "Loyalty Level": "Tier",
        # },
        # use_container_width=True,
        height=500,
    )

with col3:
    st.write("Inventory Data")
    st.dataframe(
        df3,
        # column_config={
        #     "Image": st.column_config.ImageColumn("Image", width="small"),
        #     "Price": st.column_config.NumberColumn("Price", format="$%.2f"),
        #     "Stock": st.column_config.ProgressColumn("Stock Level", min_value=0, max_value=1000),
        #     "Rating": st.column_config.NumberColumn("Rating", format="%.1f⭐"),
        #     "Stock Status": "Status",
        #     "Tags": st.column_config.ListColumn("Tags"),
        # },
        # use_container_width=True,
        height=500,
    )

# First row of charts
col1, col2, col3, col4, col5, col6 = st.beta_columns(6)

with col1:
    st.subheader("Sales Trend")
    sales_trend = df1.groupby("Quarter")["Sales"].sum().reset_index()
    # st.line_chart(sales_trend, x="Quarter", y="Sales", height=250)

with col2:
    st.subheader("Customer Segments")
    chart = alt.Chart(df2.sample(200)).mark_circle().encode(
        x=alt.X('Age', scale=alt.Scale(zero=False)),
        y=alt.Y('Total Spend', scale=alt.Scale(zero=False)),
        color='Loyalty Level',
        size='Satisfaction',
        tooltip=['Age', 'Total Spend', 'Loyalty Level', 'Satisfaction']
    ).properties(height=250)
    st.altair_chart(chart, use_container_width=True)

with col3:
    st.subheader("Stock by Category")
    stock_by_cat = df3.groupby("Category")["Stock"].sum().reset_index()
    # st.bar_chart(stock_by_cat, x="Category", y="Stock", height=250)

with col4:
    st.subheader("Sales by Quarter")
    chart_data = df1.groupby(['Quarter', 'Region'])['Sales'].mean().reset_index()
    chart = alt.Chart(chart_data).mark_circle(size=100).encode(
        x='Quarter',
        y='Sales',
        color='Region',
        size='Sales',
        tooltip=['Quarter', 'Region', 'Sales']
    ).properties(height=250)
    st.altair_chart(chart, use_container_width=True)

with col5:
    st.subheader("Price vs Rating")
    # st.scatter_chart(df3[["Price", "Rating"]].sample(100), height=250)

with col6:
    st.subheader("Growth Over Time")
    st.area_chart(df1[["Growth"]].rolling(50).mean(), height=250)

# Second row of charts
col1, col2, col3, col4, col5, col6 = st.beta_columns(6)

with col1:
    st.subheader("Price-Rating Heat")
    chart = alt.Chart(df3.sample(500)).mark_rect().encode(
        x=alt.X('Price:Q', bin=alt.Bin(maxbins=20)),
        y=alt.Y('Rating:Q', bin=alt.Bin(maxbins=20)),
        color=alt.Color('count()', scale=alt.Scale(scheme='viridis')),
        tooltip=['count()', 'Price', 'Rating']
    ).properties(height=250)
    st.altair_chart(chart, use_container_width=True)

with col2:
    st.subheader("Stock Levels")
    st.line_chart(df3["Stock"].rolling(50).mean(), height=250)

with col3:
    st.subheader("Growth Analysis")
    growth_data = df1.copy()
    growth_data['MA_Growth'] = growth_data['Growth'].rolling(50).mean()
    growth_data['MA_Profit'] = growth_data['Profit'].rolling(50).mean()
    growth_data = growth_data.dropna()

    base = alt.Chart(growth_data.sample(200)).encode(
        x=alt.X('MA_Growth:Q', title='Growth Rate (%)'),
        y=alt.Y('MA_Profit:Q', title='Profit'),
        color='Region'
    ).properties(height=250)

    points = base.mark_circle()
    lines = base.transform_regression(
        'MA_Growth', 'MA_Profit'
    ).mark_line()
    
    chart = (points + lines)
    st.altair_chart(chart, use_container_width=True)

with col4:
    st.subheader("Profit Trend")
    st.area_chart(df1["Profit"].rolling(50).mean(), height=250)

with col5:
    st.subheader("Age vs Spend")
    # st.scatter_chart(df2[["Age", "Total Spend"]].sample(100), height=250)

with col6:
    st.subheader("Regional Trends")
    region_trends = df1.groupby(['Region', 'Quarter'])['Sales'].mean().reset_index()
    region_pivot = region_trends.pivot(index='Quarter', columns='Region', values='Sales')
    # st.line_chart(
    #     region_pivot,
    #     height=250,
    #     use_container_width=True
    # )

# st.html("""
#     <style>
#     .stMainBlockContainer {
#         padding-left: 1.5rem !important; 
#         padding-right: 1.5rem !important; 
#         padding-top: 4rem !important;
#     }
#     </style>
#     """
# )